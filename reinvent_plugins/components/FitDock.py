#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import CDPL.Chem as Chem
import CDPL.Pharm as Pharm
import CDPL.ConfGen as ConfGen

import os
import sys
import shutil
import subprocess
from multiprocessing import Pool
import argparse

# =========================
# run commamd
# =========================
def run_cmd(cmd, timeout=300):
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout
        )
        if result.returncode != 0:
            print(f"[CMD ERROR]\n{cmd}\n{result.stderr.decode()}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"[TIMEOUT]\n{cmd}")
        return False


# =========================
# SMILES to SDF
# =========================
def smile2mol(smile="", num_confs=1, output_path=""):
    try:
        writer = Chem.MolecularGraphWriter(output_path)
        mol = Chem.parseSMILES(smile.strip())

        prot_state_gen = Chem.ProtonationStateStandardizer()
        ConfGen.prepareForConformerGeneration(mol)

        struct_gen = ConfGen.StructureGenerator()
        struct_gen.generate(mol)
        struct_gen.setCoordinates(mol)

        prot_state_gen.standardize(
            mol,
            Chem.ProtonationStateStandardizer.PHYSIOLOGICAL_CONDITION_STATE
        )

        Chem.perceiveComponents(mol, True)

        conf_gen = ConfGen.ConformerGenerator()
        conf_gen.settings.minRMSD = 1.0
        conf_gen.settings.energyWindow = 50
        conf_gen.settings.maxNumOutputConformers = num_confs

        conf_gen.generate(mol)
        conf_gen.setConformers(mol)

        writer.write(mol)
        writer.close()

        return True
    except Exception as e:
        print(f"smile2mol error: {e}")
        return False


# =========================
# sdf to mol2
# =========================
def sdf2mol(input="", output_path=""):
    cmd = f"obabel -isdf {input} -O {output_path}"
    return run_cmd(cmd, timeout=60)


# =========================
# mol2 to sdf
# =========================
def mol2mol(input="", output_path=""):
    try:
        reader = Chem.MoleculeReader(input)
        mol = Chem.BasicMolecule()
        reader.read(mol)

        writer = Chem.MolecularGraphWriter(input)
        writer.write(mol)
        writer.close()

        cmd1 = f"obabel -imol2 {input} -O {input} -d"
        cmd2 = f"obabel -imol2 {input} -O {output_path} --minimize --ff uff -a -p 7.2"

        ok1 = run_cmd(cmd1, 60)
        ok2 = run_cmd(cmd2, 120)

        return ok1 and ok2

    except Exception as e:
        print(f"mol2mol error: {e}")
        return False


# =========================
# FitDock
# =========================
def fitdock(ref_ligand_mol2="", query_ligand_mol2="", output_mol2=""):
    dir_ = os.path.dirname(output_mol2)
    cmd = f"cd {dir_} && FitDock -Tlig {ref_ligand_mol2} -Qlig {query_ligand_mol2} -align_only -o {output_mol2}"
    return run_cmd(cmd, timeout=300)


# =========================
# gnina minimization
# =========================
def mini_by_gnina(output_dir, input_protein_pdb, dock_sdf_path, ref_ligand_sdf):

    energy_min_mol2 = os.path.join(output_dir, "postdock_min.mol2")
    energy_min_log = os.path.join(output_dir, "postdock_min.log")

    cmd = f"gnina1.3.1 -r {input_protein_pdb} -l {dock_sdf_path} --autobox_ligand {ref_ligand_sdf} --log {energy_min_log} -o {energy_min_mol2} --minimize --cnn_scoring none"# --local_only 
    
    if not run_cmd(cmd, timeout=600):
        return None, 0.0, 10.0

    energy_min_sdf = os.path.join(output_dir, "postdock_min.sdf")
    mol2mol(input=energy_min_mol2, output_path=energy_min_sdf)

    score = 0.0
    RMSD_ = 10.0

    if os.path.exists(energy_min_log):
        with open(energy_min_log, "r") as f:
            for line in f:
                if line.startswith("Affinity"):
                    score = float(line.split()[1])
                elif line.startswith("RMSD"):
                    RMSD_ = float(line.split()[1])

    return energy_min_sdf, score, RMSD_


# =========================
# main
# =========================
def fitdock_single(input_smi, input_protein_pdb, ref_ligand_mol2, output_dir, ddG, rmsd):

    try:
        os.makedirs(output_dir, exist_ok=True)

        query_sdf = os.path.join(output_dir, "query.sdf")
        if not smile2mol(input_smi, output_path=query_sdf):
            return 1.0

        query_mol2 = os.path.join(output_dir, "query.mol2")
        if not sdf2mol(query_sdf, query_mol2):
            return 1.0

        dock_mol2 = os.path.join(output_dir, "dock.mol2")

        rec_copy = os.path.join(output_dir, "protein.pdb")
        shutil.copy(input_protein_pdb, rec_copy)

        ref_copy = os.path.join(output_dir, "ref_lig.mol2")
        shutil.copy(ref_ligand_mol2, ref_copy)

        if not fitdock(ref_copy, query_mol2, dock_mol2):
            return 1.0

        dock_sdf = os.path.join(output_dir, "dock.sdf")
        if not mol2mol(dock_mol2, dock_sdf):
            return 1.0

        ref_sdf = os.path.join(output_dir, "ref_lig.sdf")
        if not os.path.exists(ref_sdf):
            mol2mol(ref_copy, ref_sdf)

        _, refscore, _ = mini_by_gnina(output_dir, rec_copy, ref_sdf, ref_sdf)
        _, score, RMSD_ = mini_by_gnina(output_dir, rec_copy, dock_sdf, ref_sdf)

        if RMSD_ < rmsd:
            result = score - refscore
        else:
            result = 1.0

        return result

    except Exception as e:
        print(f"[MAIN ERROR] {e}")
        return 1.0


# =========================
# worker
# =========================

def fitscore(
        output_dir,
        smiles,
        protein_pdb,
        ref_ligand_mol2,
        ddg,
        rmsd
):
    os.makedirs(output_dir, exist_ok=True)

    fitdock_input = []
    for i, smi in enumerate(smiles):
        if smi.strip():
            fitdock_input.append(
                (smi.strip(),
                    protein_pdb,
                    ref_ligand_mol2,
                    os.path.join(output_dir, str(i)),
                    ddg,
                    rmsd)
            )

    # multiprocess
    process_num=max(1, os.cpu_count() - 5)
    with Pool(processes=process_num, maxtasksperchild=10) as pool:
        results = pool.starmap(fitdock_single, fitdock_input)
    
    return results
    

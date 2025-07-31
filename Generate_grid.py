#!/usr/bin/env python3
"""
Automated Synthetic Spectrum Generation Pipeline
===============================================
Generates synthetic stellar spectra using SPECTRUM (v2.77c) by Richard Gray
over a grid of stellar parameters: Teff, logg, [M/H], vsini

Author: Generated for stellar parameter extraction pipeline
"""

import os
import sys
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path

# Constants
c = 299792458.0  # Speed of light in m/s (for rotational broadening)

# Parameter grid definition - More realistic ranges
TEFF_RANGE = np.arange(4000, 7001, 250)
LOGG_RANGE = np.arange(1.0, 5.1, 0.5)
MH_VALUES = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.2, 0.5]
VSINI_VALUES = [0, 5, 10, 20, 50]

user_input = input("Enter start, stop and step values for wavelength in commas: ")
WAVE_START, WAVE_END, WAVE_STEP = map(float, user_input.split(','))

SUPERMODEL_MAP = {
    -2.5: "am25k2odfnew.dat",
    -2.0: "am20k2odfnew.dat",
    -1.5: "am15k2odfnew.dat",
    -1.0: "am10k2odfnew.dat",
    -0.5: "am05k2odfnew.dat",
     0.0: "ap00k2odfnew.dat",
     0.2: "ap02k2odfnew.dat",
     0.5: "ap05k2odfnew.dat"
}

def check_spectrum_installation():
    print("\U0001F527 Checking SPECTRUM installation...")
    executables = ['selectmod', 'spectrum', 'avsini']
    for exe in executables:
        if os.path.exists(exe) and os.access(exe, os.X_OK):
            print(f"✅ Found {exe} in current directory")
            continue
        try:
            subprocess.run(['which', exe], check=True, capture_output=True)
            print(f"✅ Found {exe} in PATH")
        except subprocess.CalledProcessError:
            print(f"❌ Error: {exe} not found!")
            return False
    return True

def find_executable(exe_name):
    if os.path.exists(exe_name) and os.access(exe_name, os.X_OK):
        return f"./{exe_name}"
    try:
        result = subprocess.run(['which', exe_name], check=True, capture_output=True, text=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        pass
    common_paths = [
        f"/Users/utkarsh8022/Downloads/spectrum277c/{exe_name}",
        f"./spectrum277c/{exe_name}"
    ]
    for path in common_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            return path
    raise FileNotFoundError(f"Could not find executable: {exe_name}")

def setup_directories():
    for d in ['extracted_models', 'Synfiles']:
        Path(d).mkdir(exist_ok=True)
    print("✓ Output directories ready")

def check_input_files():
    print("\U0001F527 Checking input files...")
    if not os.path.exists('luke.lst'):
        print("❌ Error: line list file 'luke.lst' not found")
        return False
    if not os.path.exists('Supermodels'):
        print("❌ Error: 'Supermodels' directory not found")
        return False
    missing_models = []
    for mh, filename in SUPERMODEL_MAP.items():
        filepath = f"Supermodels/{filename}"
        if not os.path.exists(filepath):
            missing_models.append(f"[M/H]={mh}: {filename}")
    if missing_models:
        print("❌ Missing supermodel files:")
        for model in missing_models:
            print(f"   {model}")
        return False
    return True

def estimate_microturbulence(teff, logg):
    return round(1.1 + 1.6e-4 * (teff - 5500) - 0.3 * (logg - 4.0), 2)

def get_supermodel_path(mh):
    if mh not in SUPERMODEL_MAP:
        raise ValueError(f"No supermodel available for [M/H] = {mh}")
    model_path = f"Supermodels/{SUPERMODEL_MAP[mh]}"
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Supermodel not found: {model_path}")
    return model_path

def extract_model_atmosphere(supermodel_path, teff, logg, mh):
    mh_str = f"p{abs(mh):02.0f}" if mh >= 0 else f"m{abs(mh)*10:02.0f}"
    output_mod = f"extracted_models/{mh_str}t{teff:04.0f}g{logg:.1f}.mod"
    selectmod_exe = find_executable('selectmod')
    cmd = [selectmod_exe, supermodel_path, output_mod, str(teff), str(logg)]
    subprocess.run(cmd, check=True)
    if not os.path.exists(output_mod):
        raise FileNotFoundError(f"selectmod did not create {output_mod}")
    return output_mod

def generate_spectrum(mod_file, vturb, teff, logg, mh):
    spectrum_exe = find_executable('spectrum')
    mh_str = f"p{abs(mh):02.0f}" if mh >= 0 else f"m{abs(mh)*10:02.0f}"
    output_path = f"Synfiles/{mh_str}t{teff:04.0f}g{logg:.1f}.spc.r00"

    spectrum_input = f"""{mod_file}
luke.lst
{output_path}
{vturb}
{WAVE_START},{WAVE_END}
{WAVE_STEP}
"""

    print(f"   ⏳ Running SPECTRUM for Teff={teff}, logg={logg}, [M/H]={mh}")
    print(f"   → Output file: {output_path}")

    process = subprocess.Popen(
        [spectrum_exe, 'n'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    try:
        stdout, stderr = process.communicate(input=spectrum_input, timeout=300)

        if process.returncode != 0:
            raise RuntimeError(f"SPECTRUM failed with code {process.returncode}: {stderr}")

        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise RuntimeError(f"No output spectrum created: {output_path}")

        return output_path

    except subprocess.TimeoutExpired:
        process.kill()
        raise RuntimeError("SPECTRUM execution timed out")

def apply_rotational_broadening(base_file, teff, logg, mh):
    avsini_exe = find_executable('avsini')
    mh_str = f"p{abs(mh):02.0f}" if mh >= 0 else f"m{abs(mh)*10:02.0f}"
    broadened_files = []
    for vsini in VSINI_VALUES:
        if vsini == 0:
            broadened_files.append(base_file)
            continue
        output_file = base_file.replace(".spc.r00", f".spc.r{vsini:02d}")
        cmd = [avsini_exe, base_file, output_file, str(vsini), "0.6", str(WAVE_STEP)]
        subprocess.run(cmd, check=True)
        broadened_files.append(output_file)
    return broadened_files

def main():
    print("\n🌟 Starting SPECTRUM synthetic grid generation...")
    if not check_spectrum_installation():
        return 1
    if not check_input_files():
        return 1
    setup_directories()

    combinations = [(teff, logg, mh) for teff in TEFF_RANGE for logg in LOGG_RANGE for mh in MH_VALUES]
    total_spectra = len(combinations) * len(VSINI_VALUES)

    print(f"\n📊 Total combinations: {len(combinations)}")
    print(f"📁 Total spectra to be generated (including vsini): {total_spectra}")
    confirm = input("Proceed with generation? [Y/n]: ").strip().lower()
    if confirm not in ['', 'y', 'yes']:
        print("Aborted by user.")
        return 0

    for teff, logg, mh in tqdm(combinations, desc="Generating spectra"):
        try:
            print(f"\n⚙ Generating: Teff={teff} K, logg={logg}, [M/H]={mh:+.1f}")
            mod_file = extract_model_atmosphere(get_supermodel_path(mh), teff, logg, mh)
            vturb = estimate_microturbulence(teff, logg)
            base_file = generate_spectrum(mod_file, vturb, teff, logg, mh)
            apply_rotational_broadening(base_file, teff, logg, mh)
        except Exception as e:
            print(f"❌ Failed for Teff={teff}, logg={logg}, [M/H]={mh} → {e}")

if __name__ == "__main__":
    sys.exit(main())


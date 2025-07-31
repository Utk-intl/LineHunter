import subprocess
from pathlib import Path
import shutil
import re


base_dir = Path("/Users/utkarsh8022/Downloads/spectrum277c")
spectrum_exec = base_dir / "spectrum"
line_list = "luke.lst"
mod_source_dir = Path("/Users/utkarsh8022/Downloads/marcs_st_krz/marcs_atm_models")
output_spc_dir = base_dir / "marcs_synfiles"
wavelength_range = (3000.0, 4800.0)
wavelength_step = 0.01
timeout_limit = 300  # seconds


def is_model_safe(teff, logg, feh, vt):
    return (
        teff >= 3000 and
        0.5 <= logg <= 5.5 and
        feh >= -4.0 and
        vt >= 0.0
    )


def parse_params(filename):
    try:
        teff = float(re.search(r"([ps])(\d{4})", filename).group(2))
        logg = float(re.search(r"g([+-]\d+\.\d+)", filename).group(1))
        vt_str = re.search(r"t(\d\d)", filename).group(1)
        vt = float(vt_str[0] + "." + vt_str[1])
        feh = float(re.search(r"z([+-]\d+\.\d+)", filename).group(1))
        return teff, logg, feh, vt
    except Exception as e:
        print(f"⚠️ Could not parse parameters from {filename}: {e}")
        return None


print("📂 Copying .mod files to spectrum277c ...")
mod_files = sorted(mod_source_dir.glob("*.mod"))
if not mod_files:
    print("⛔ No .mod files found in source directory.")
    exit(1)

for mod_file in mod_files:
    shutil.copy(mod_file, base_dir / mod_file.name)
print(f"✅ Copied {len(mod_files)} files to {base_dir.name}/")


output_spc_dir.mkdir(parents=True, exist_ok=True)

atm_files = sorted(base_dir.glob("*.mod"))

for atm_path in atm_files:
    fname = atm_path.name
    params = parse_params(fname)
    if not params:
        continue

    teff, logg, feh, vt = params

    if not is_model_safe(teff, logg, feh, vt):
        print(f"⛔ Skipping: {fname} → Unsafe parameters Teff={teff}, logg={logg}, [M/H]={feh}, vt={vt}")
        continue

    output_name = fname.replace(".mod", ".spc")
    output_path = output_spc_dir / output_name

    print(f"✅ Processing {fname} → {output_name} | Teff={teff}, logg={logg}, vt={vt}")

    try:
        proc = subprocess.run(
            [str(spectrum_exec)],
            input=f"{fname}\n{line_list}\n{output_name}\n{vt}\n{wavelength_range[0]},{wavelength_range[1]}\n{wavelength_step}\n",
            capture_output=True,
            text=True,
            cwd=base_dir,
            timeout=timeout_limit
        )

        if "now exiting to system" in proc.stdout or "series failed in expint" in proc.stdout:
            print(f"❌ Error running SPECTRUM on {fname}\n↪ {proc.stdout.splitlines()[-3:]}")
        elif not output_path.exists():
            fallback_path = base_dir / output_name
            if fallback_path.exists():
                shutil.move(fallback_path, output_path)
                print(f"📁 Moved {output_name} → marcs_synfiles/")
            else:
                print(f"❌ Output .spc file not found for {fname}")
        elif output_path.stat().st_size == 0:
            print(f"❌ Empty .spc file generated for {fname} ❗ Check input model or vt.")
        else:
            print(f"🎉 Done: {output_name} ({output_path.stat().st_size / 1024:.1f} KB)")

    except subprocess.TimeoutExpired:
        print(f"🔥 Timeout: SPECTRUM took too long on {fname}")
    except Exception as e:
        print(f"💥 Unexpected error on {fname}: {e}")


print("🧹 Cleaning up copied .mod files ...")
for mod_file in base_dir.glob("*.mod"):
    try:
        mod_file.unlink()
        print(f"✅ Deleted {mod_file.name}")
    except Exception as e:
        print(f"⚠️ Could not delete {mod_file.name}: {e}")
print("🧹 Cleanup complete.")


spc_files = list(output_spc_dir.glob("*.spc"))
if not spc_files:
    print("⛔ No .spc files were generated. Please check the input models and parameters.")
    exit(1)
else:
    print(f"✅ {len(spc_files)} .spc files generated in {output_spc_dir.name}/")
    print("🎯 All tasks completed successfully!")

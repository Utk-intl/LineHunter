import os
kerg = 1.38054e-16

def read_marcs_krz(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        teff_line = lines[1]
        teff = float(teff_line.split('EFF=')[1].split()[0])
        logg = float(teff_line.split('GRAV=')[1].split()[0])

        mh = 0.0
        vturb = 0.0
        for part in filename.split('_'):
            if part.startswith('z'):
                mh = float(part[1:])
            if part.startswith('t'):
                try:
                    vturb = float(part[1:])
                except ValueError:
                    vturb = 0.0

        layers = None
        for line in lines:
            parts = line.strip().split()
            if parts and parts[-1].isdigit():
                val = int(parts[-1])
                if 10 < val < 200:
                    layers = val
                    break

        if layers is None:
            raise ValueError("Couldn't determine number of layers.")

        print(f'Teff : {teff}, Logg : {logg}, [M/H] : {mh}, Layers : {layers}')
        return lines, teff, logg, mh, vturb, layers

def extract_structure_data(lines, layers,vturb):
    data = []
    count = 0
    kerg = 1.38054e-16
    for line in lines:
        if ',' in line:
            parts = [x.strip() for x in line.split(',')]
            if len(parts) >= 5:
                try:
                    rhox = float(parts[0])
                    temp = float(parts[1])
                    ne = float(parts[2])
                    ntmne = float(parts[3])
                    rho = float(parts[4])
                    pg = (ntmne+ne)*kerg*temp
                    junk = 0.0
                    vt = vturb * 1.0e+5
                    data.append([rhox, temp, pg, ne, junk, junk, vt])
                    count += 1
                    if count >= layers:
                        break
                except ValueError:
                    continue
    if len(data) != layers:
        raise ValueError("Incorrect number of layers parsed.")
    return data

# NEW: Batch loop over all .krz files
krz_files = sorted([f for f in os.listdir('.') if f.endswith('.krz')])
output_dir = input("Enter the name of the folder you want to save the .mod files in (default='marcs_atm_models')>").strip()
if not output_dir:
    output_dir = 'marcs_atm_models'

os.makedirs(output_dir, exist_ok=True)
print(f'✅ All .mod files will be saved in: {output_dir}')
if not krz_files:
    print("No .krz files found in the current directory.")
else:
    for model in krz_files:
        try:
            print(f"\n📄 Processing: {model}")
            model_base = os.path.splitext(model)[0]
            model_output = os.path.join(output_dir, model_base + '.mod')

            lines, teff, log, mh, vturb, layers = read_marcs_krz(model)
            atm_data = extract_structure_data(lines, layers, vturb)

            g = open(model_output, "w")
            g.write("%8.1f %f %f %d\n" % (teff, log, mh, layers))
            for i in range(layers):
                g.write("%11.9e %8.1f %7.5e %7.5e %7.5e %7.5e %7.5e\n" %
                        (atm_data[i][0], atm_data[i][1], atm_data[i][2],
                         atm_data[i][3], atm_data[i][4], atm_data[i][5], atm_data[i][6]))
            g.close()
            print(f"✅ Saved: {model_output}")

        except Exception as e:
            print(f"❌ Error processing {model}: {e}")


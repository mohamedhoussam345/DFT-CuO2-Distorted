import numpy as np
import os

# CONSTANTS
# ==============================
hbar_au   = 1.0
Ry_to_eV  = 13.605693
cm1_to_au = 4.55634e-6    # 1 cm-1 in atomic units
amu_to_me = 1822.888486

# SYSTEM PARAMETERS
# ==============================
alat  = 7.2188            # in bohr
Omega = 590.8464          # unit cell volume in a.u.^3
Ncell = 1                 # primitive cell
grid  = np.array([45, 45, 216])
Ngrid = np.prod(grid)
dV    = Omega / Ngrid

# Atomic masses in a.u.
M = np.array([63.546, 16.0, 16.0]) * amu_to_me  # Cu, O, O

# Atomic positions in alat units (from ATOMIC_POSITIONS)
R_tau = np.array([
    [0.0302820, 0.0302820, 0.0000000],  # Cu
    [0.5209351, 0.0410225, 0.0000000],  # O
    [0.0410225, 0.5209351, 0.0000000],  # O
])

base = r"D:\final project"

# TRIPLETS: (k_idx, kp_idx, q_idx)
# ==============================
K_LIST = [1, 20, 27, 29, 30]
triplets = [(ik, ikp, 1) for ik in K_LIST for ikp in K_LIST]

# q-points in crystal coords 
# ==============================
q_cryst = np.array([
    [ 0.0000,  0.0000, 0.0],  # q1
    [ 0.0000,  0.1250, 0.0],  # q2
    [ 0.0000,  0.2500, 0.0],  # q3
    [ 0.0000,  0.3750, 0.0],  # q4
    [ 0.0000, -0.5000, 0.0],  # q5
    [ 0.1250,  0.1250, 0.0],  # q6
    [ 0.1250,  0.2500, 0.0],  # q7
    [ 0.1250,  0.3750, 0.0],  # q8
    [ 0.2500, -0.5000, 0.0],  # q9
    [ 0.2500, -0.3750, 0.0],  # q10
    [ 0.2500, -0.2500, 0.0],  # q11
    [ 0.2500, -0.1250, 0.0],  # q12
    [ 0.2500,  0.2500, 0.0],  # q13
    [ 0.2500,  0.3750, 0.0],  # q14
    [ 0.3750, -0.5000, 0.0],  # q15
    [ 0.3750, -0.3750, 0.0],  # q16
    [ 0.3750, -0.2500, 0.0],  # q17
    [ 0.3750,  0.3750, 0.0],  # q18
    [-0.5000, -0.5000, 0.0],  # q19
    [-0.5000, -0.3750, 0.0],  # q20
    [-0.3750, -0.3750, 0.0],  # q21
])

# LOAD WAVEFUNCTION
# ==============================
WFC_FILES = {
    1:  "kpoint1.txt",
    20: "kpoint20.txt",
    27: "kpoint27.txt",
    29: "kpoint29.txt",
    30: "kpoint30.txt",
}

def load_psi(kidx):
    path = os.path.join(base, WFC_FILES[kidx])
    psi = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("(") and line.endswith(")"):
                line = line[1:-1]
                parts = line.split(",")
                psi.append(float(parts[0]) + 1j*float(parts[1]))
    return np.array(psi)

# LOAD DVSCF MODE (txt, two columns: real imag)
# skip first 3 modes (acoustic), skip negative freq modes
# ==============================
def load_dvscf_mode(iq, imode):
    path = os.path.join(base, f"dvscf_q{iq}_mode{imode}.txt")
    data = np.loadtxt(path, comments='#')
    return data[:, 0] + 1j * data[:, 1]  # flatten complex array

# LOAD POLARIZATION VECTORS from dynmat_7q.npy
# ===============================
dynmat = np.load(os.path.join(base, "dynmat_7q.npy"), allow_pickle=True).item()

def load_polarization(iq):
    """Returns freqs (nmodes,) and eigvecs list of (nat,3) complex arrays."""
    data = dynmat[iq]
    return data['freqs'], data['eigvecs']

# BUILD r-GRID
# =============================
ix, iy, iz = np.meshgrid(
    np.arange(grid[0]),
    np.arange(grid[1]),
    np.arange(grid[2]),
    indexing="ij"
)
r_frac = np.stack([
    ix.flatten() / grid[0],
    iy.flatten() / grid[1],
    iz.flatten() / grid[2]
], axis=1)  # shape (Ngrid, 3)

#  wavefunctions
psi_cache = {}
def get_psi(kidx):
    if kidx not in psi_cache:
        psi_cache[kidx] = load_psi(kidx)
    return psi_cache[kidx]

# MAIN LOOP
# ============================
SKIP_MODES = 3  # skip first 3 acoustic modes

results = []
print(f"{'k':>4} {'k_prime':>7} {'q':>4} {'mode':>5} {'freq(cm-1)':>12} "
      f"{'|D| (eV)':>12} {'Re(D)':>14} {'Im(D)':>14}")
print("-" * 90)

for (ik, ikp, iq) in triplets:
    psi_k  = get_psi(ik)
    psi_kp = get_psi(ikp)

    # q vector in crystal coords
    q_vec = q_cryst[iq - 1]

    # load polarization for this q
    freqs_cm1, eigvecs = load_polarization(iq)

    for imode in range(9):
        # skip first 3 acoustic modes
        if imode < SKIP_MODES:
            continue

        freq = freqs_cm1[imode]

        # skip imaginary frequencies
        if freq <= 0:
            print(f"  Skipping k{ik},k'{ikp},q{iq} mode {imode+1}: "
                  f"freq={freq:.3f} cm-1 (imaginary)")
            continue

        omega_au = freq * cm1_to_au

        # load dvscf for this mode (already mode-projected)
        dV_mode = load_dvscf_mode(iq, imode + 1)  # 1-based file naming

        # polarization vector for this mode: shape (nat=3, 3)
        e_nu = eigvecs[imode]  # (3, 3)

        # sum over atoms tau
        D = 0.0 + 0.0j
        for tau in range(3):
            # mass weighted prefactor for this atom
            prefactor = np.sqrt(hbar_au / (2 * M[tau] * Ncell * omega_au))

            # norm of polarization vector for this atom
            e_norm = np.linalg.norm(e_nu[tau])

            # phase factor: exp(i q . R_tau) in crystal coords
            phase_tau = np.exp(1j * 2 * np.pi * np.dot(q_vec, R_tau[tau]))

            # real-space integral: <psi_k'| dV_mode |psi_k>
            integrand = np.conj(psi_kp) * dV_mode * psi_k
            integral  = (dV / Omega) * np.sum(integrand)

            D += prefactor * e_norm * phase_tau * integral

        D *= Ry_to_eV  # convert to eV

        results.append({
            'k': ik, 'kp': ikp, 'q': iq,
            'mode': imode + 1, 'freq': freq, 'D': D
        })

        print(f"{ik:>4} {ikp:>7} {iq:>4} {imode+1:>5} {freq:>12.3f} "
              f"{np.abs(D):>12.6e} {D.real:>14.6e} {D.imag:>14.6e}")

# SAVE RESULTS
# ==============================
out_txt = os.path.join(base, "g_matrix_elements.txt")
with open(out_txt, 'w', encoding='utf-8') as f:
    f.write(f"{'k':>4} {'k_prime':>7} {'q':>4} {'mode':>5} {'freq(cm-1)':>12} "
            f"{'|D|(eV)':>12} {'Re(D)':>14} {'Im(D)':>14}\n")
    f.write("-" * 80 + "\n")
    for r in results:
        f.write(f"{r['k']:>4} {r['kp']:>7} {r['q']:>4} {r['mode']:>5} "
                f"{r['freq']:>12.3f} {np.abs(r['D']):>12.6e} "
                f"{r['D'].real:>14.6e} {r['D'].imag:>14.6e}\n")

print(f"\nTotal matrix elements computed: {len(results)}")
print(f"Saved to {out_txt}")

# NORMALIZATION CHECK
# ==============================
psi_test = get_psi(1)
norm = np.sum(np.abs(psi_test)**2) * dV / Omega
print(f"\nNormalization check |psi_k1|^2 dV/Omega = {norm:.6f}  (should be ~1.0)")


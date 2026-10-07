# 2D Poisson Solver with HDF5 Storage and ParaView Export

Finite-difference solution of the 2D Poisson equation on the unit square, with a convergence study,
results stored in a hierarchical **HDF5** file and exported to **VTK** (`.vtu`) for visualisation in **ParaView**.

Master 1 CHPS (High-Performance Computing and Simulation), Université de Perpignan Via Domitia,
course *Outils pour la Visualisation* (2025–2026).

| Approximate solution (ParaView, N = 10) | Convergence study |
|---|---|
| ![Solution for N = 10 displayed in ParaView](figures/paraview_N10.png) | ![Error and time versus mesh size](figures/convergence.png) |

## The problem

$$-\Delta u = f \quad \text{on } (0,1)^2, \qquad u = 0 \text{ on the boundary}, \qquad f(x,y) = \sin(2\pi x)\sin(2\pi y)$$

The exact solution is $u(x,y) = \dfrac{\sin(2\pi x)\sin(2\pi y)}{8\pi^2}$.

The domain is discretised with $N \times N$ interior nodes ($h = 1/(N+1)$) and the 5-point stencil.
The sparse matrix is assembled with Kronecker products and solved with SciPy's sparse direct solver.
The error is the discrete L2 norm $h\,\lVert U - U_{\text{exact}} \rVert_2$.

## Results

| N | h = 1/(N+1) | L2 error | Time (s) |
|---:|---:|---:|---:|
| 10 | 9.09e-02 | 1.75e-04 | 0.003 |
| 20 | 4.76e-02 | 4.75e-05 | 0.002 |
| 50 | 1.96e-02 | 8.02e-06 | 0.008 |
| 100 | 9.90e-03 | 2.04e-06 | 0.033 |

The observed order of convergence (log-log regression) is **2.007**, as expected for a second-order scheme.
Times were measured on the machine that generated the files in `results/`; they will vary on yours.

**Performance note.** The first version of this code factorised the *dense* matrix (`lu_factor(A.toarray())`).
Switching to a sparse solver, and replacing the Python loops by vectorised assembly, brings
N = 100 from about 18 s to about 0.03 s on the same machine, with identical results.

## Quick start

```bash
git clone https://github.com/<your-username>/poisson2d-hdf5-vtk.git
cd poisson2d-hdf5-vtk
pip install -r requirements.txt

python run_convergence.py      # solve, write results/resultats.hdf5, make the figures
python export_vtu.py           # HDF5 -> results/solution_N10.vtu, open it in ParaView
```

Options:

```bash
python run_convergence.py --sizes 10 20 50 100 200 400
python export_vtu.py --group Sol_N50 --out results/solution_N50.vtu
```

In ParaView: open the `.vtu` file, click *Apply*, colour by `U`, and use *Surface With Edges* to see the mesh.

## HDF5 layout

```text
resultats.hdf5
├── Sol_N10/            attribute: N = 10
│   └── Solution        dataset: U at the N*N interior nodes (numbered k = i + j*N)
├── Sol_N20/ ...
├── Sol_N50/ ...
├── Sol_N100/ ...
├── N_values            dataset: [10, 20, 50, 100]
├── L2_error            dataset: L2 error for each N
└── CPU_Time            dataset: time in seconds for each N
```

Inspect it with `h5dump results/resultats.hdf5` or `h5py`.

## Repository layout

```text
poisson2d.py          solver: assembly, solve, exact solution, convergence order
run_convergence.py    runs several N, writes the HDF5 file and the figures
export_vtu.py         reads the HDF5 file and writes a VTK unstructured grid (.vtu)
results/              HDF5 file and a sample .vtu file (N = 10)
figures/              convergence plot, solution plot, ParaView screenshot
```

## Author

Komlan Gagnon. MIT license, see [LICENSE](LICENSE).

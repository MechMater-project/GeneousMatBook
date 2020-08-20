# User Read-me

The GeneousMatbook is a set of jupyterlab notebooks designed to calculate the homogeneous effective behaviors of heterogeneous materials :
- Matrix with spherical and ellipsoidal isotropic and anisotropic fillers

It is designed to:
- calculate effective mechanical properties, 
- Provide comprehensive descriptions of the implemented models,
- give an easy environment to introduce new models

## Installation
Some sections of the code require you to have gfortran installed on your computer and available on your Path. Depending on your environnement, here are the steps you should follow when running the software for the fisrt time.

### Issue with gfortran and windows Path : 
Mac/windows: You may find a good gfortran compiler on this website https://gcc.gnu.org/wiki/GFortranBinaries

Linux: apt install gfortran

If you have issue with adding gfortran with your Path on windows, this website may help https://davescience.wordpress.com/2012/03/21/gfortran-on-windows/. **Furthermore you should install mingw64 in C: problems occurs when your binaries are in a file with space in the name or too far from root.**

### Local with pip:
- Install JupyterLab
- Open a terminal on JupyterLab
- Exectute the following commands
    - pip3 install ipywidgets
    - pip3 install nodejs
    - pip3 install ipympl
    - jupyter labextension install @jupyter-widgets/jupyterlab-manager
    - jupyter labextension install jupyter-matplotlib
    - jupyter nbextension enable --py widgetsnbextension
- Restart JupyterLab

### Local with conda
- Install Anaconda
- Open a new terminal on JupyterLab
- Execute the following commands
    - conda install -c conda-forge ipywidgets
    - conda install -c conda-forge ipympl
    - conda install -c conda-forge nodejs
    - jupyter labextension install @jupyter-widgets/jupyterlab-manager
    - jupyter lab build
- Restart Anaconda and JupyterLab

### Compiling locally Fortran functions with f2py
- Required the first time 'homogenization_main.ipynb' is used
- Open a new terminal on JupyterLab, and run:
    - f2py -c fortran_tools.f -m fortran_tools
- For developpers, this step will be performed when changing fortran_tools.f, which may never happen.

## Using the software

The software enables you to compute the homogenized behaviors of simple microstructures using models implemented with Python. The two notebooks 'homogenization_main.ipynb' and 'homogenization_visco.ipynb' contain a user-friendly friendly interface. On the former, the user will find multiple sections:

- Computation of the homogenized behavior of defined microstructures

The section first invites the user to manually build inclusions by selecting their geometry and their behavior. Note that the inclusion is not generated as long as the 'Generate inclusion' button is not toggled. Once the inclusions generated, the user can create a microstructure by adding and removing any of the previously generated inclusions. The microstructure is generated after the 'Generate microstructure' button is toggled. The script will then display the compatible models and enable the user to chose one for the computation of the homogenized behavior.

The model comparison enables the user to compare multiple models in terms of properties with respect to the volume fraction of fillers for the last generated microstructure. 


- Automatic computation from a .txt file

Enables the computation of multiple homogenized behaviors from a text file. The input file must have a .txt extension and placed in the inputs/automatic_computation folder. example.txt gives an exemple for the expected input format. Its first line must be '*homogenization'. Two different microstructure must be separated with a line containing the '*' sign.

*homogenization

Self-consistent <- first line is the name of the model

K:15,G:15       <- second line stands for the istropic behavior of the matrix it can be 'K:value,G:value' for bulk modulus and shear modulus or 'E:value,nu:value' for Young modulus and Poisson's ratio

0,1.0,0.3,K:200,G:200 <- third line describes the first inclusion, start with '0,1.0' for spheres, then the volume fraction of filler and its mechanical behavior. For ellipsoid fillers the two first values '0,1.0' are changed into its aspect ratios 'a2/a1,a3/a1'.
For an anisotropic inclusion 'C' indicates that the user is given the stiffness tensor and 'S' the compliance tensor and the name of the file where C or S is stored is given. For instance,

10,100,0.2,C,cristalPET_Cij.txt 

'*   ' <- line to end the microstructure.  Note that several inclusions may be entered before this line when considering a matrix filled with different inclusions.

In the case of the 4-phase model where the inclusion is surrounded by an interphase, the inclusion + interphase data are given as

*4-phase self-consistent

K:15,G:15

0.,0.2,K:200,G:200,0.3,K:1000,G:10 <- where '0.2,K:200,G:200' are the volume fraction and the behavior of the filler  and '0.3,K:1000,G:10' are the volume fraction and the behavior of the interphase


- Model description

Displays the descriptions of the implemented models. The descriptions are Markdown files in the inputs/model_comparison folder. If the Latex equations are not properly displayed (known bug), the user should click on the blue bar on the left of the description twice and the latex equations will be computed.

- Estimates of parameters by an inverse method

Computes optimal microstructure parameters to reach a given target homogenized behavior. If the optimization problem admits multiple local minima, the program will only return one of them (not necessarily the best). Runs well for a few unknown parameters (volume fraction of an inclusion, or for a single unknown matrix behavior parameter), but can still be improved.

## Differences :


- 'homogenization_main.ipynb' and 'homogenization_visco.ipynb'

The first handle only elastic cases, whereas the second can compute linear viscoelastic homogenized behavior for **single spherical isotropic inclusion**.

The homogenization_visco.ipynb will be released soon
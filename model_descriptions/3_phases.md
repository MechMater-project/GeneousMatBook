# 3-phase or generalized self-consistent (GSC) model 
[R.M., Christensen, K.H., Lo, 1979. Solution for effective shear properties in three phase sphere and cylinder models. J. Mech. Phys. Solids, 27, 315-330. Erratum, 1986. J. Mech. Phys. Solids, 34, 639]

Self-consistent solution for Hashin composite-spheres randomly dispersed

Notation:
$G$ shear modulus, 
$K$ Bulk modulus, 
$\nu$ Poisson's ratio, 
$E$ Young modulus, 
volume fraction of fillers $f$,
index $f$ for filler, index $m$ for matrix, and $\bar{X}$ for the homogeneous medium.
## Spherical isotropic inclusions ramdomly dispersed
\begin{equation}
\bar{K}=K_m+\frac{f}{\frac{1}{K_f-K_m}+\frac{3(1-f)}{3K_m+4G_m}}\,.
\end{equation}
The shear modulus $\bar{G}$ of the composite is the positive root of:
\begin{equation}
A\,\bar{G}^2-B\,G_m\,\bar{G}+C\,G_m^2=0
\end{equation}

where
\begin{gather*}
A=8(5\nu_m -4)ge_1\,f_p^{10/3}+D
+50(8\nu_m ^2-12\nu_m +7)ge_2\,f_p+4(10\nu_m -7)e_2\,e_3\,,\\
B=4(5\nu_m -1)ge_1\,f_p^{10/3}+2D
-150(\nu_m -3)\nu_m g\,e_2\,f_p+3(15\nu_m -7)e_2\,e_3\,,\\
C=-4(5\nu_m -7)ge_1\,f_p^{10/3}+D
-25(\nu_m ^2-7)ge_2\,f_p+(5\nu_m +7)e_2\,e_3
\end{gather*}
with
\begin{gather*}
g=G_f/G_m-1\,,\quad
D=2\left[63ge_2+2e_1\,e_3\right]f_p^{7/3}-252ge_2\,f_p^{5/3}\,,\\
e_1=(49+35\nu_f-70\nu_m -50\nu_f\,\nu_m )g+105(\nu_f-\nu_m )\,,\\
e_2=(7+5\nu_f)g+35(1-\nu_f)\,,\quad
e_3=2(4-5\nu_m )g+15(1-\nu_m )\,.
\end{gather*}

## Model validation
This model has been validated on data  Fig 3 p12 in [Christensen R.M., and Lo K.H.,  1979. SOLUTIONS FOR EFFECTIVE SHEAR PROPERTIES IN THREE PHASE SPHERE AND CYLINDER MODELS],

<img src="model_descriptions/model_validate/GSC_Christensen_G.png" alt="drawing" width="500">

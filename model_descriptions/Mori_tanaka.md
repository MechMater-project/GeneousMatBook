# Mori-Tanaka
[Mori, T. and K. Tanaka, 1973. Average Stress in Matrix and Average Elastic Energy of Materials with Misfitting Inclusions. Acta Metallurgica, 21(5), 571-574.]

Notation:
$G$ shear modulus, 
$K$ Bulk modulus, 
$\nu$ Poisson's ratio, 
$E$ Young modulus, $C$ elastic stiffness tensor,
volume fraction of fillers $f_f$, volume fraction of matrix $f_m$
index $f$ for filler, index $m$ for matrix, and $\bar{X}$ for the homogeneous medium.
## General equations

$$\bar{C}=C_m+f_f(C_f-C_m):A^{MT}$$

with $C_m$ the matrix behavior $C_f$ the filler behavior, $f$ its volume fraction, and finally $A^{MT}=A^{Esh}\left( f_m\,I +f_f A^{Esh}\right)^{-1}$

with $A^{Esh}=\left( I+S(C_m^{-1}C_f-I) \right)^{-1}$ the Eshelby localization tensor defined by
the Eshelby tensor $S(C_m)$ that depends on the inclusion geometry and the matrix behaviour only.

## Spherical isotropic inclusions ramdomly dispersed

Inclusions are assumed surrounding by a layer of matrix but some of the inclusion interactions are taken into account.

$$\bar{K}=K_m+\frac{f(K_f-K_m)(3K_m+4G_m)}{3K_f+4G_m+3(1-f)(K_f-K_m)}$$
$$\bar{G}=G_m+\frac{5fG_m(G_f-G_m)(3K_m+4G_m)}{5G_m(3K_m+4G_m)+6(1-f)(G_f-G_m)(K_m+2G_m)}$$

## Ellipsoidal inclusions ramdomly dispersed
We consider the general case of anisotropic inclusions randomly dispersed in an isotropic matrix. The resulting equivalent homogeneous material is isotropic. 
The Eshelby tensor $S$ expressions are given in [Mura, T. Micromechanics of defects in solids, second revised edition, Dordecht: Kluwer, 1987.]

$n$ inclusions are randomly dispersed in space, and each inclusion $i$ is considered as a filler constitutive phase of volume fraction $f_i=f_f/n$. The consitutive equations writes now as:
$$\bar{C}=C_m+\sum_i^n f_i(C_{f_i}-C_m):A_i^{MT}$$ 
with $C_{f_i}$ and $A_i^{MT}=\left( f_m\,I +f_i A^{Esh}\right)^{-1}$ the behavior of the filler and of the self-consistent localization tensor computed in the material reference. 

## Model validation
This model have been validated on data Eq. (64) p.15 in [Christensen,  R.M.,  1990.  A critical evaluation for a class of micromechanics models], 

<img src="model_descriptions/model_validate/MT_Christensen_G1.png" alt="drawing" width="500">
<img src="model_descriptions/model_validate/MT_Christensen_K1.png" alt="drawing" width="500">
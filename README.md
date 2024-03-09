# Exciton in WSe2 quantum dot using configuration interaction

## Bethe Salpeter equation for an exciton

Configuration interaction (CI) algorithm can be used to study the problem of an exciton when it is formulated in the basis of all possible excitations 

$|km>=c_m^{\dagger}c_k|GS>=c_m^{\dagger}c_k \prod_{p filled} c_{p\uparrow}^{\dagger}c_{p\downarrow}^{\dagger}|0>$,

where $p,q,k,m$ denote single particle (SP) states. The single exciton wavefunction then reads

$|X>=\sum_{km\sigma}A_{km\sigma}|km\sigma>$,

where the coefficients $A_{km\sigma}$ are obtained by soving the Bethe Salpeter equation 

$(\epsilon_{m\sigma}+\Sigma_{m\sigma}-\epsilon_{k\sigma}-\Sigma_{k\sigma})A_{km\sigma}+\sum_{ij\sigma'}(V_{kjim}-V_{kjmi})A_{ij\sigma'}=EA_{km\sigma}$

in the matrix form, which can be built by a CI algorithm in the basis of excitations. $\Sigma$ is self energy and $V_{kjim}$ etc. are scattering Coulomb matrix elements (CME).


## Input files

The input files consist of
* SP energy levels for spin up (U) and down (D)
* CME elements for spin combinations UU, DD, UD
* quantum numbers specific to the problem of MoS2 QD: valley index and angular momentum

## References

[1] D. Miravet, L.Szulakowska,  M. Bieniek, M. Korkusinski, P. Hawrylak, ”Theory of excitonic complexes in gated WSe2 quantum dots”, (in preparation) 2024




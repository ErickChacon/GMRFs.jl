# Gaussian Markov random fields

Implementation of GMRF for spatial analysis.


## Regular grids

![regular](figures/01-regular-grid.png)

## Irregular latices

![irregular](figures/02-irregular-grid.png)

## Graphs

![graphs](figures/03-graph.png)

## Todo

- [ ] `GCMRF`:
    - [ ] Reparametrize the `struct` using the `base` of the precision matrix.
    - [ ] Finish implementation of `logpdf!` using `DFT` and `IDFT`.

## Compute the derivatives of Gram-Schmidt vectors

Gram-Schmidt orthogonalization is a well-known technique where an non-orthogonal set of column vectors is made orthogonal via

![equation](https://latex.codecogs.com/svg.latex?%5Cbegin%7Balign*%7D%20%7B%5Cbf%20w%7D_i%20%3D%20%7B%5Cbf%20q%7D_i%20-%20%5Csum_%7Bj%3D1%7D%5E%7Bi-1%7D%5Cleft%28%5Cfrac%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20q%7D_i%7D%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20w%7D_j%7D%5Cright%29%7B%5Cbf%20w%7D_j%2C%20%5Cquad%20i%20%3D%201%2C%5Ccdots%2C%20d.%20%5Cend%7Balign*%7D)

This repository provides a Python3 script for computing the exact derivatives ![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cpartial%7B%5Cbf%20w%7D_i/%5Cpartial%7B%5Cbf%20q%7D_k)  and  ![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cpartial%28%7B%5Cbf%20w%7D_i/%5CVert%7B%5Cbf%20w%7D_i%5CrVert_2%29/%5Cpartial%7B%5Cbf%20q%7D_k).

The following files are present:

* `gram_schmidt_derivatives.py`: a Python3 class for computing Gram-Schmidt vectors and the derivatives.
* `derivation.pdf`: a document which details the derivation of the derivative expressions.
* `check_derivatives.ipynb`: a Jupyter notebook which verifies the validity of the derivative expressions, and which verifies the Python3 implementation of `gram_schmidt_derivatives.py`. Also acts as a tutorial for the Python3 class.

Development was funded by the EU Horizon 2020 project [VECMA](http://www.vecma.eu/).

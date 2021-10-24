## Differentiating Gram Schmidt

Gram-Schmidt orthogonalization is a well-known technique where an non-orthogonal set of column vectors ![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%7B%5Cbf%20q%7D_i%5Cin%5Cmathbb%7BR%7D%5ED) is made orthogonal via

![equation](https://latex.codecogs.com/svg.latex?%5Cbegin%7Balign*%7D%20%7B%5Cbf%20w%7D_i%20%3D%20%7B%5Cbf%20q%7D_i%20-%20%5Csum_%7Bj%3D1%7D%5E%7Bi-1%7D%5Cleft%28%5Cfrac%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20q%7D_i%7D%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20w%7D_j%7D%5Cright%29%7B%5Cbf%20w%7D_j%2C%20%5Cquad%20i%20%3D%201%2C%5Ccdots%2C%20d.%20%5Cend%7Balign*%7D)

We have developed a simple recurrence relation to compute the derivatives of ![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%7B%5Cbf%20w%7D_i) with respect to ![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%7B%5Cbf%20q%7D_i):

![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cboxed%7B%5Cbegin%7Balign*%7D%20%5Cfrac%7B%5Cpartial%7B%5Cbf%20w%7D_1%7D%7B%5Cpartial%7B%5Cbf%20q%7D_1%7D%20%3D%3A%20D_%7B11%7D%20%3D%20I_D%20%5C%5C%20%5Cfrac%7B%5Cpartial%7B%5Cbf%20w%7D_i%7D%7B%5Cpartial%7B%5Cbf%20q%7D_i%7D%20%3D%3A%20D_%7Bii%7D%20%3D%20D_%7Bi-1%2C%5C%2C%20i-1%7D%20-%20%5Cfrac%7B%7B%5Cbf%20w%7D_%7Bi-1%7D%7B%5Cbf%20w%7D_%7Bi-1%7D%5ET%7D%7B%7B%5Cbf%20w%7D_%7Bi-1%7D%5ET%7B%5Cbf%20w%7D_%7Bi-1%7D%7D%2C%5Cquad%20i%20%3E%201%2C%20%5C%5C%20%5Cfrac%7B%5Cpartial%7B%5Cbf%20w%7D_i%7D%7B%5Cpartial%7B%5Cbf%20q%7D_k%7D%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7Bi-1%7DD_%7Bij%7D%5Cfrac%7B%5Cpartial%7B%5Cbf%20w%7D_j%7D%7B%5Cpartial%7B%5Cbf%20q%7D_k%7D%2C%5Cquad%20i%5Cneq%20k%2C%20%5Cquad%20i%20%3E%20k%2C%20%5C%5C%20D_%7Bij%7D%3A%3D%20-%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%7B%5Cbf%20w%7D_j%7D%5Cleft%5B%5Cleft%28%5Cfrac%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20q%7D_i%7D%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20w%7D_j%7D%5Cright%29%7B%5Cbf%20w%7D_j%5Cright%5D%20%3D%20-%5Cleft%5B%5Cfrac%7B1%7D%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20w%7D_j%7D%5C%3B%7B%5Cbf%20w%7D_j%7B%5Cbf%20q%7D_i%5ET%20-%20%5Cfrac%7B2%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20q%7D_i%7D%7B%28%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20w%7D_j%29%5E2%7D%5C%3B%7B%5Cbf%20w%7D_j%7B%5Cbf%20w%7D_j%5ET%20&plus;%20%5Cfrac%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20q%7D_i%7D%7B%7B%5Cbf%20w%7D_j%5ET%7B%5Cbf%20w%7D_j%7D%5C%3B%20I_D%5Cright%5D%2C%5Cquad%20i%5Cneq%20j%2C%5Cquad%20i%20%3E%20j.%20%5Cend%7Balign*%7D%7D)

Gram-Schmidt vectors are often normalized in practise, to obtain the set:

![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cbegin%7Balign*%7D%20W%20%3D%20%5Cbigg%5C%7B%5Cfrac%7B%7B%5Cbf%20w%7D_1%28%7B%5Cbf%20q%7D_1%29%7D%7B%5ClVert%7B%5Cbf%20w%7D_1%28%7B%5Cbf%20q%7D_1%29%5CrVert_2%7D%2C%5C%3B%5C%3B%5Cfrac%7B%7B%5Cbf%20w%7D_2%28%7B%5Cbf%20q%7D_1%2C%20%7B%5Cbf%20q%7D_2%29%7D%7B%5ClVert%7B%5Cbf%20w%7D_2%28%7B%5Cbf%20q%7D_1%2C%20%7B%5Cbf%20q%7D_2%29%5CrVert_2%7D%2C%5C%3B%5C%3B%5Ccdots%5C%3B%5C%3B%2C%5Cfrac%7B%7B%5Cbf%20w%7D_d%28%7B%5Cbf%20q%7D_1%2C%7B%5Cbf%20q%7D_2%5C%2C%5Ccdots%2C%7B%5Cbf%20q%7D_d%29%7D%7B%5ClVert%7B%5Cbf%20w%7D_d%28%7B%5Cbf%20q%7D_1%2C%7B%5Cbf%20q%7D_2%5C%2C%5Ccdots%2C%7B%5Cbf%20q%7D_d%29%5CrVert_2%7D%5Cbigg%5C%7D.%20%5Cend%7Balign*%7D)

The derivatives of these normalized vectors are obtained by premultiplying the unnormalized derivative with a matrix that only depends upon ![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%7B%5Cbf%20w%7D_i):

![equation](https://latex.codecogs.com/svg.latex?%5Cinline%20%5Cbegin%7Balign*%7D%20%5Cboxed%7B%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%7B%5Cbf%20q%7D_k%7D%5Cleft%28%5Cfrac%7B%7B%5Cbf%20w%7D_i%7D%7B%5ClVert%7B%5Cbf%20w%7D_i%5CrVert_2%7D%5Cright%29%20%3D%20%5Cleft%5B%5Cfrac%7BI_D%7D%7B%5ClVert%7B%5Cbf%20w%7D_i%5CrVert_2%7D%20-%20%5Cfrac%7B%7B%5Cbf%20w%7D_i%7B%5Cbf%20w%7D_i%5ET%7D%7B%5ClVert%7B%5Cbf%20w%7D_i%5CrVert%5E3_2%7D%5Cright%5D%5Cfrac%7B%5Cpartial%7B%5Cbf%20w%7D_i%7D%7B%5Cpartial%7B%5Cbf%20q%7D_k%7D.%7D%20%5Cend%7Balign*%7D)

### Contents repository

The following files are present in the repository:

* `gram_schmidt_derivatives.py`: a Python3 class for computing Gram-Schmidt vectors and the derivatives.
* `derivation.pdf`: a document which details the derivation of the derivative expressions.
* `check_derivatives.ipynb`: a Jupyter notebook which verifies the validity of the derivative expressions, and which verifies the Python3 implementation of `gram_schmidt_derivatives.py`. Also acts as a tutorial for the Python3 class.

Development was funded by the EU Horizon 2020 project [VECMA](http://www.vecma.eu/).

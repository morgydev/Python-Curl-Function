# Python-Curl-Function
An extension to numpy using discrete fourier transforms to compute the curl of 2D and 3D functions. This produces results far more accurate than using 10th-order finite difference derivatives (*which is as far as I tested*).

## Main Functions
```curl2D``` computes the curl for a 2D input function $\boldsymbol{F} = (F_x, F_y)$.

```curl3D``` computes the curl for a 3D input function $\boldsymbol{F} = (F_x, F_y, F_z)$.

```rcurl2D``` computes the curl for a 2D purely real input function $\boldsymbol{F} = (F_x, F_y)$.

```rcurl3D``` computes the curl for a 3D input real function $\boldsymbol{F} = (F_x, F_y, F_z)$.

## Usage
![alt text](https://imgur.com/a/AnMvL3l)

## Explanation
Mathematically, the curl $\nabla \times \boldsymbol{F}$ can be written as

$`\begin{vmatrix}
     \hat{i} & \hat{j} & \hat{k}\\ 
     \frac{\partial}{\partial x} & \frac{\partial}{\partial y} & \frac{\partial}{\partial z}\\
     F_x(x, y, z) & F_y(x, y, z) & F_z(x, y, z)
\end{vmatrix}`$

which can be further expanded to

$= \left( \frac{\partial F_z}{\partial y} - \frac{\partial F_y}{\partial z} \right)\hat{i} + \left( \frac{\partial F_x}{\partial z} - \frac{\partial F_z}{\partial x} \right)\hat{j} + \left( \frac{\partial F_y}{\partial x} - \frac{\partial F_x}{\partial y} \right)\hat{k}$

Clearly, there are 6 derivatives (2 in the 2D case) required to compute. This code utilises the discrete fourier transform to compute these, which (based on testing) are more accurate than at least 10th-order finite difference derivatives, at the cost of being more computationally expensive.

$\hat{f}(\vec{k}) = \frac{1}{N} \Sigma_{m = 0}^{N-1} f(\vec{x}) e^{-i\vec{k} \bullet \vec{x}_m}$

$f(\vec{x}) = \Sigma_{k = 0}^{N-1} \hat{f}(\vec{k}) e^{i\vec{k} \bullet \vec{x}_m}$

where $\vec{x}_m = -2 \pi \vec{x}/N$

## Further Optimisation
For real functions, you can instead use ```rcurl2D``` and ```rcurl3D``` to perform a real discrete fourier transform, which is Hermitian-symmetric and thus we can discard half of the transformed values to save on memory and computation.

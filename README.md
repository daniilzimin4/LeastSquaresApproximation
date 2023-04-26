## LeastSquaresApproximation

# Jacobi:
α = I - (D_1 * A)\
β = D_1 * b\
X0 = β\
X(i+1) = X(i) + D_1 * (b - A*X(i) )

such that:\
D : the diagonal matrix which contains the diagonal elements from A\
D_1 : the inverse matrix of D\
I : the identity matrix\
A, b : they are given

# Seidel:
α = I - (D_1 * A)\
β = D_1 * b\
X0 = β\
X(i+1) = B_1 * (b - C*X(i))\
A, b : they are given

such that:
D : the diagonal matrix which contains the diagonal elements from A\
D_1 : the inverse matrix of D\
B : the lower triangular part of A (with the diagonal)\
C : the upper triangular part of A (without the diagonal)\
B_1 : the inverse matrix of B

# Epsilon
epsilon could be calculated by finding the norm of the vector X(i+1) - X(i)

# Gravity

An object at position $\mathbf r_0$ will experience a vector acceleration $\mathbf a$ due to another body with mass $m_1$ at $\mathbf r_1$

$$
\mathbf a = \frac{G m_1}{|\mathbf r_1- \mathbf r_0|^3}(\mathbf r_1 - \mathbf r_0)
$$
where 

If there are multiple bodies, we add the accelerations vectorially:
$$
\mathbf a  = \sum_{i=1}^N \frac{G m_i}{|\mathbf r_i- \mathbf r_0|^3}(\mathbf r_i- \mathbf r_0)
$$

Python code that implements different approximation methods using object oriented techniques for an arbitrary number of bodies, given an initial position, velocity and a mass for each body.
The appromiation methods used are; Euler's Method, Runge-Kutta Method 2 and Runge-Kutta Method 4

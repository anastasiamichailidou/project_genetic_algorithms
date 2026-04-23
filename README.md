Project: Genetic Algorithms for Function Approximation

DescriptionThis project was developed for the Optimization Techniques course (7th Semester, December 2025).
The objective is to find an analytical expression for an unknown continuous function of two variables $f(u_1, u_2)$.To approximate the function, a linear combination of up to 15 Gaussian functions is used.
The optimization of the Gaussian parameters (centers $c$ and deviations $\sigma$) is performed using a Genetic Algorithm implemented from scratch in Matlab, without using built-in libraries like ga().

Technical SpecificationsModel: Linear combination of Gaussian functions.
Target Function: $f(u_1, u_2) = \sin(u_1 + u_2) \sin(u_2^2)$ for evaluation purposes.
Domain: $u_1 \in [-1, 2]$ and $u_2 \in [-2, 1]$.
Evaluation: Tested on a separate set of input-output data from the one used for identification.

Project Structure

Matlab_Files/: Contains all the source code for the project.
Report.pdf: Technical report including problem description, results, and observations.

Requirements
Matlab (R2025 version or later).

Submission
The project was completed for the February 2026 examination period.

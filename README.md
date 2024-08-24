# Optical Ray Tracer Simulation

## Project Overview

This project involves the development of a simulation to analyze and optimize optical lens configurations, with a primary focus on minimizing spherical aberration. The simulation compares the performance of various lens types, including plano-convex and optimized lenses, by tracing light rays through these optical systems and evaluating their impact on image formation.

## Key Features

- **Ray and Optical Elements Classes:**
  - Developed classes to simulate the propagation of light rays through different optical elements, including spherical refracting surfaces.
  - Implemented methods for calculating ray-surface intercepts and refracted directions using Snell's law in vector form.

- **Ray Tracing and Visualization:**
  - Created functions to trace rays through the optical system and visualize their paths, including options to display spot diagrams and calculate the RMS (Root Mean Square) radius at the paraxial focus.

- **Lens Optimization:**
  - Developed an optimization routine using `scipy.optimize` to adjust lens curvatures for minimizing the RMS radius, thereby reducing spherical aberration.
  - Compared the performance of optimized lenses with standard plano-convex lenses as well as with a commercial lens from ThorLabs.

## Results

### **Spherical Aberration Reduction:**
- **Plano-Convex vs. Optimized Lens:**
  - The optimized lens configuration, with a combination of positive and smaller negative curvatures, demonstrated a significant reduction in spherical aberration compared to the plano-convex lens.
  - The RMS radius for the optimized lens was reduced by up to 75% in certain configurations, resulting in a more precise focus and improved image quality.

### **Comparison with ThorLabs Lens:**
- **Simulation of ThorLabs LBF254-075-A Lens:**
  - The simulation was extended to replicate the ThorLabs LBF254-075-A lens, which is designed to minimize spherical aberration with specific radii of curvature.
  - The optimized lens produced by the simulation showed similar radii of curvature and performance characteristics to the ThorLabs lens.
  
- **Performance Metrics:**
  - Although the ThorLabs lens slightly outperformed the optimized lens in terms of RMS radius (0.0127 mm vs. 0.0129 mm), the results suggest that the optimization function successfully approached the commercial standard.
  - The minor performance gap may be attributed to the optimization process reaching a local minimum rather than a global minimum, indicating potential areas for further refinement.

## Conclusion

This project successfully developed a robust simulation tool for analyzing and optimizing lens configurations in optical systems. The findings demonstrate that custom lens optimizations can approach the performance of commercial lenses like those from ThorLabs, with potential for further improvement through enhanced optimization techniques.

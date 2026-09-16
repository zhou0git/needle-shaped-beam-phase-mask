# Needle-Shaped Beam Phase Pattern Generator

This MATLAB program generates an SLM phase pattern for a needle-shaped beam. It spatially multiplexes the phases of multiple axial foci, iteratively adjusts their positions to improve axial intensity uniformity, and adds a linear phase to control lateral displacement in the focal plane.

## Main Features

- Configure the optical system, focal points, and SLM parameters.
- Quasi-randomly assign SLM pixels to different axial foci.
- Generate and iteratively optimize the needle-shaped beam phase.
- Add a two-dimensional linear phase for lateral beam displacement.
- Embed the active phase region into the full SLM frame.
- Display and save the resulting grayscale BMP phase pattern.
- Report the final focus positions and relative center intensities.

## Usage

Run the following script in MATLAB:

```matlab
needle_beam_phase_mask
```

Optical parameters, the number of foci, the phase adjustment coefficient, the iteration count, and SLM settings can be modified in the `Parameter settings` section.

The displacement magnitude and direction can be changed in `Add a linear phase`. The output directory and filename can be changed in `Generate, display, and save the SLM phase pattern`.

## Code Structure

- `Parameter settings`: Defines the parameters and creates the SLM coordinate grid.
- `Needle-shaped beam phase`: Generates and optimizes the multiplexed phase.
- `Add a linear phase`: Controls lateral displacement in the focal plane.
- `Generate, display, and save the SLM phase pattern`: Quantizes, embeds, displays, and saves the phase pattern.
- `Iterative algorithm`: Updates the focus positions and calculates their center intensities.
- `Pixel allocation function`: Randomly assigns unit-cell pixels to the foci.

## Output

The program generates a full-size 8-bit grayscale SLM phase pattern and saves it as a BMP file. The MATLAB workspace also retains the continuous phase, quantized phase, focus positions, and center-intensity results.

## Notes

- The number of foci must be a perfect square.
- Pixel allocation is randomized, so results may vary slightly between runs. Set a fixed random seed when reproducibility is required.
- Ensure that the output directory exists or that MATLAB has permission to create it.
- Large parameter settings may require significant memory for the three-dimensional pixel-allocation array.

## Attribution

This repository contains an independent MATLAB implementation of the needle-shaped beam generation method described in the cited publication. It is not the original implementation released by the authors of that publication.

The underlying method and its scientific concepts are attributed to the original authors. Users of this software should cite the original publication when this implementation contributes to academic or scientific work.

## Disclaimer

This repository is an independent research implementation based on the cited publication. It is not affiliated with, authorized by, or endorsed by the original authors or their institutions.

The software is provided for research and educational purposes. Users are responsible for verifying the implementation, parameters, numerical results, and suitability for their own optical systems and applications.

The authors and contributors of this repository make no guarantee regarding the accuracy, completeness, reliability, or fitness of the software for any particular purpose and accept no liability for any loss or damage arising from its use.

## Reference

This implementation is based on the method described in:

J. Zhao, Y. Winetraub, L. Du, A. Van Vleck, K. Ichimura, C. Huang, S. Z. Aasi, K. Y. Sarin, and A. de la Zerda, “Flexible method for generating needle-shaped beams and its application in optical coherence tomography,” *Optica*, vol. 9, no. 8, pp. 859–867, 2022. [https://doi.org/10.1364/OPTICA.456894](https://doi.org/10.1364/OPTICA.456894)

If this code contributes to published work, please cite the original article above.

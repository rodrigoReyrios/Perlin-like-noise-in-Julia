Perlin-like Noise in Julia
=============

<!-- ![banner]() -->
<!-- ![badge]() -->
<!-- ![badge]() -->
Simple implementation of a Perlin noise algorithm as described in [these](https://www.cs.umd.edu/class/spring2018/cmsc425/index.shtml) course notes by Dave Mount and Roger Eastman at the University of Maryland. These notes are in turn based on Ken Perlin's improved algorithm archived [here](https://web.archive.org/web/20180418002912/http://mrl.nyu.edu/~perlin/paper445.pdf).

Contents
-----

`noise.jl` exports the following:
- Function to generate pseudo-random vectors across a 2D patch

- Function to generate Perlin noise across a 2D mesh

![perlin mesh visual](/figs/PerlinNoiseDemo.png)

See `noise_illustration.ipynb` for an illustration of the algorithm and a small demo on the how to generate the above image.

Author
-------
* Rodrigo Rios [:email:](rodrigoreyrios@gmail.com)
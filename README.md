# NLA

## Description

**Disclaimer:** This is very much a work in progress. Use at your own risk and please report any problems to Leonardo Capitani **(**leocapi07\@gmail.com**).**

This repository contains data and code to implement Nested Library Analysis (NLA) in R by using [Gaussian Process-EDM](https://tanyalrogers.github.io/GPEDM) method (Munch, S. & Rogers, T. 2025). NLA was proposed by [Huang et al. 2024](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011759#sec009): *Detecting shifts in nonlinear dynamics using Empirical Dynamic Modeling with Nested-Library Analysis*.

### Background: 

NLA was presented to answer this research question: Does the variable of interest present any change point in time?

[Huang et al. 2024](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011759#sec009) found that NLA is capable to detect change points and outperforms existing approaches based on statistical characteristics. Indeed, [Huang et al. 2024](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011759#sec009) work shows that shifts in mean and variance cannot serve as a generic signal of regime shifts.

NLA repository associated to the original publication is stored [here](https://github.com/chromian/Nested-Library-Analysis).

### Goal

We aim to translate the NLA core algorithms into R and we want to actualize the Empirical Dynamic Modelling prediction method from S-Map to Gaussian Processes (GP).

### License

MIT License © 2026 Leonardo Capitani

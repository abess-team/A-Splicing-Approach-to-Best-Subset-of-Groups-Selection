
# Group splicing numerical experiments
This repository contains scripts to run the synthetic datasets and real-world dataset analysis described
in *A Splicing Approach to Best Subset of Groups Selection*. 

## Codes

* `Synthetic_dataset_analysis/Synthetic_dataset_analysis.R` : R script used to run the synthetic datasets analysis.
* `Real-world_dataset_analysis/Real-world_dataset_analysis.R` : R script used to run the real-world dataset analysis.
* `Real-world_dataset_analysis/trim32.rda` : Dataset used in the real-world dataset analysis. 
* `gomp/gomp.R` : R interface of `gomp.cpp`.
* `gomp/gomp.cpp` : C++ implementation of group orthogonal matching pursuit (GOMP).

## Softwares

* Group Lasso : R package `grpreg` (3.4.0).
* Group MCP : R package `grpreg` (3.4.0).
* GOMP : Implementation in R language with `Rcpp` modules.
* Group Splicing : R package `abess` (0.4.0).


## Citations

Please cite the following publications if you make use of the material here.

- Yanhang Zhang, Junxian Zhu, Jin Zhu, and Xueqin Wang (2022). A Splicing Approach to Best Subset of Groups Selection. arXiv preprint arXiv:2104.12576.

- Jin Zhu, Liyuan Hu, Junhao Huang, Kangkang Jiang, Yanhang Zhang, Shiyun Lin, Junxian Zhu, and Xueqin Wang (2021). abess: A Fast Best Subset Selection Library in Python and R. arXiv preprint arXiv:2110.09697, 2021.

The corresponding BibteX entries:

```
@article{2104.12576,
author = {Yanhang Zhang and Junxian Zhu and Jin Zhu and Xueqin Wang},
title = {A Splicing Approach to Best Subset of Groups Selection},
year = {2022},
journal = {arXiv preprint arXiv:2104.12576},
}
```
and
```
@article{zhu-abess-arxiv,
  author    = {Jin Zhu and Liyuan Hu and Junhao Huang and Kangkang Jiang and Yanhang Zhang and Shiyun Lin and Junxian Zhu and Xueqin Wang},
  title     = {abess: A Fast Best Subset Selection Library in Python and R},
  journal   = {arXiv preprint arXiv:2110.09697},
  year      = {2021},
}
```


## Contact
Please direct questions and comments to the [issues page](https://github.com/abess-team/Group-splicing_codes/issues).

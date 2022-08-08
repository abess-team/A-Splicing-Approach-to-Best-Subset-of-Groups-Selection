
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

- Yanhang Zhang, Junxian Zhu, Jin Zhu, and Xueqin Wang (2021). Certifiably Polynomial Algorithm for Best Group Subset Selection. arXiv preprint arXiv:2104.12576, 2021.

- Jin Zhu, Xueqin Wang, Liyuan Hu, Junhao Huang, Kangkang Jiang, Yanhang Zhang, Shiyun Lin and Junxian Zhu (2022). abess: A Fast Best-Subset Selection Library in Python and R. Journal of Machine Learning Research, 23(202), 1-7.

The corresponding BibteX entries:

```
@article{zhang2021certifiably,
  title = {Certifiably polynomial algorithm for best group subset selection},
  author = {Yanhang Zhang and Junxian Zhu and Jin Zhu and Xueqin Wang},
  journal = {arXiv preprint arXiv:2104.12576},
  year = {2021}
}
```
and
```
@article{JMLR:v23:21-1060,
  author  = {Jin Zhu and Xueqin Wang and Liyuan Hu and Junhao Huang and Kangkang Jiang and Yanhang Zhang and Shiyun Lin and Junxian Zhu},
  title   = {abess: A Fast Best-Subset Selection Library in Python and R},
  journal = {Journal of Machine Learning Research},
  year    = {2022},
  volume  = {23},
  number  = {202},
  pages   = {1--7},
  url     = {http://jmlr.org/papers/v23/21-1060.html}
}
```


## Contact
Please direct questions and comments to the [issues page](https://github.com/abess-team/Group-splicing_codes/issues).

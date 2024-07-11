Download the data collection we used from the [website](https://zenodo.org/records/10642766)
The downloaded `crossdocked_pocket10.zip` is ready for unzipping in the `data` folder. The unzipped
data in the folder should be like this format:

```bash
DATA
└─crossdocked_pocket10
    ├─1433B_HUMAN_1_240_pep_0
    ├─1433C_TOBAC_1_256_0
    ├─...
└─crossdocked_pocket10_name2id.pt
└─crossdocked_pocket10_processed.lmdb
└─split_by_name.pt
└─test_list.tsv
```

With the respect to the source of dataset, if this data is helpful, please also cite:
```bash
@article{francoeur2020three,
  title={Three-dimensional convolutional neural networks and a cross-docked data set for structure-based drug design},
  author={Francoeur, Paul G and Masuda, Tomohide and Sunseri, Jocelyn and Jia, Andrew and Iovanisci, Richard B and Snyder, Ian and Koes, David R},
  journal={Journal of chemical information and modeling},
  volume={60},
  number={9},
  pages={4200--4215},
  year={2020},
  publisher={ACS Publications}
}
```
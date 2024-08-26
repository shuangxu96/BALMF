# Hyperspectral Image Denoising by Asymmetric Noise Modeling
This repository provides the code and implementation details for the paper "Hyperspectral Image Denoising by Asymmetric Noise Modeling" published in IEEE TGRS. The paper proposes a novel denoising framework, bandwise asymmetric Laplacian noise modeling matrix factorization (BALMF), that leverages the asymmetric nature of hyperspectral image noise to achieve superior denoising performance.

## Key Features
* **Asymmetric Noise Modeling**:  The framework explicitly models the asymmetric noise characteristics
* **Non-i.i.d. Noise Modeling**:  Each band is assigned a specific noise model, leading to characterize non-i.i.d. noise along spectral dimension.

## Usage
To run the code, please follow the instructions in the `README.md` file. The repository includes the following scripts:
* `BALMF.m`: BALMF model.
* `demo_GF5.m`: Evaluate the denoising performance on a GF5 dataset.

## Citation
If you use this code in your research, please cite the corresponding paper:
```
@article{BALMF,
  author       = {Shuang Xu and
                  Xiangyong Cao and
                  Jiangjun Peng and
                  Qiao Ke and
                  Cong Ma and
                  Deyu Meng},
  title        = {Hyperspectral Image Denoising by Asymmetric Noise Modeling},
  journal      = {{IEEE} Trans. Geosci. Remote. Sens.},
  volume       = {60},
  pages        = {1--14},
  year         = {2022},
  doi          = {10.1109/TGRS.2022.3227735},
}
```

## Contact
If you have any questions or need further assistance, please contact Shuang Xu at xs@nwpu.edu.cn



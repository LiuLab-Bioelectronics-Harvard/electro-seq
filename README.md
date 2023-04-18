<a name="readme-top"></a>

<h3 align="center">In situ electro-seq</h3>

  <p align="center">
    Multimodal charting of molecular and functional cell states via in situ electro-seq
    <br />
    <a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4173435"><strong>Explore the manuscript (A new version will be uploaded soon.)
</strong></a>
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#Data">Data</a>
      <ul>
      </ul>
    </li>
    <li><a href="#Citation">Citation</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

Paired mapping of single-cell gene expression and electrophysiology is essential to understand gene-to-function relationships in electrogenic tissues. Here, we developed in situ electro-sequencing (electro-seq) that combines flexible bioelectronics with in situ RNA sequencing to stably map millisecond-timescale electrical activity and profile single-cell gene expression from the same cells across intact biological networks including cardiac and neuron patches. When applied to human-induced pluripotent stem cell-derived cardiomyocyte patches, in situ electro-seq enabled multimodal in situ analysis of cardiomyocyte electrophysiology and gene expression at the cellular level, jointly defining cell states and developmental trajectories. Using machine learning-based cross-modal analysis, in situ electro-seq identified gene-to-electrophysiology relationships throughout cardiomyocyte development and accurately reconstructed the evolution of gene expression profiles based on long-term stable electrical measurements. In situ electro-seq could be applicable to create spatiotemporal multimodal maps in electrogenic tissues, potentiating the discovery of cell types and gene programs responsible for electrophysiological function and dysfunction.

### Built With

* python 3.7
* jupyter notebook



<!-- Data -->
## Data
For all cardiac spatial transcriptomics data, please go to [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1346)
Other data have been attached in this GitHub


<!-- Citation -->
## Citation

Please consider cite us if you find the data/code is useful to you.

<!-- LICENSE -->
## License

Distributed under the GPL-3.0 license. See `LICENSE.txt` for more information.





<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
We note that the data analyses for in situ electro-seq relies on many computational methods. We listed the most critical tools used in our paper here:
* For more details about spatial transcriptomics cell segmentation, please refer to [ClusterMap](https://github.com/wanglab-broad/ClusterMap)
* For more details about single-cell clustering analysis and cell typing, please refer to [Scanpy](https://scanpy.readthedocs.io/en/stable/)
* For more details about the ephys-gene joint integration analysis, please refer to [Weighted Nearest Neighbor](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html)
* For more details about the ephys-gene relationship analysis, please refer to [Sparse Reduced Rank Regression](https://github.com/berenslab/patch-seq-rrr)
* For more details about the cross-modal inference, please refer to [coupled AutoEncoder](https://github.com/AllenInstitute/coupledAE-patchseq)



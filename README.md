# rQR_t
Official code of rQR_t

We noticed a significant technical discrepancy in the nomenclature used in the post-discussion justification. The reviewer claims that 'Standard CholeskyQR3' inherently uses a regularized Gram matrix. However, according to the seminal work by Fukaya et al. (2020), this is the definition of Shifted CholeskyQR3 (SCholeskyQR3), a distinct variant designed to address cases where the standard CholeskyQR3 (which we use as our baseline) fails.

The confusion between these two distinct algorithms is a critical factual error. It leads to an incorrect assessment of our paper's originality, as we are being compared against a 'shifted' baseline that was not explicitly identified as such.

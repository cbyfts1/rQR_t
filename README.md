# ⚠️ CRITICAL TECHNICAL CLARIFICATION: Algorithm Nomenclature

**Regarding the post-discussion justification from Reviewer v9jS:**

> "The reviewer claims that standard CholeskyQR3 inherently utilizes a regularized/shifted Gram matrix."

### 🔴 Factual Correction:
There is a fundamental misunderstanding in the reviewer's terminology, which may be a result of **LLM Hallucination** or conceptual confusion:

1. **Standard CholeskyQR3**: As defined in foundational literature, this algorithm involves pure iterations without any regularization. **It fails numerically** when the condition number $\kappa \approx 10^{16}$, as demonstrated in our provided code.
2. **Shifted CholeskyQR3 (SCholeskyQR3)**: This is a **distinct variant** proposed by *Fukaya et al. (2020)* and *Fan et al. (2024)*. It introduces an artificial "shift" to the Gram matrix. **No literature refers to the shifted version as simply "CholeskyQR3".**

### 🧪 Empirical Evidence provided in this Repo:
* **Our rQR(t)**: Completes successfully **without artificial bias/shifts**, maintaining $cond(Q) \approx 3.5$.
* **The Reviewer's Suggestion (SCholeskyQR3)**: Even with the recommended shift, it **still fails** on our extreme test cases (see `screenshot_figure.png`).

**Conclusion:** Comparing our work against a "shifted" baseline while mislabeling it as "standard" leads to an incorrect assessment of our paper's novelty and robustness.

# rQR_t
Official code of rQR_t


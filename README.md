# pleiotest

## Testing of Simulated Pleiotropic Genotypes and Phenotypes

### Overview

The `pleiotest` package provides a comprehensive framework for evaluating pleiotropy in simulated genetic data. It includes multiple statistical methods for testing the association of genetic variants with multiple correlated traits, leveraging a range of model-based and data-driven approaches. 

Key features include:
- Multiple pleiotropy tests, such as `MixFisher`, `MixVar`, `MixTippett`, `MixSD`, `PCMinP`, and `PCLC`
- Integration with [`pleiosim`](https://github.com/Broccolito/pleiosim) for simulating pleiotropic genotypes and phenotypes
- Utilization of summary statistics from GWAS-like studies
- Support for correlated traits and network-wide association studies

This package is designed for researchers studying the genetic basis of complex traits and diseases.

## Installation

To install the development version of `pleiotest` from GitHub, you will need **`devtools`**:

```r
# Install devtools if not installed
install.packages("devtools")

# Install pleiotest from GitHub
devtools::install_github("Broccolito/pleiotest")
```

### Dependencies

`pleiotest` depends on several R packages, including:
- [`MPAT`](https://github.com/Broccolito/MPATclone)
- [`MTARclone`](https://github.com/Broccolito/MTARclone)
- [`pleiosim`](https://github.com/Broccolito/pleiosim)

If these are not automatically installed, you can install them manually:

```r
devtools::install_github("Broccolito/MPATclone")
devtools::install_github("Broccolito/MTARclone")
devtools::install_github("Broccolito/pleiosim")
```

## Usage

### Load the Package

```r
library(pleiotest)
```

### Example: Running a Pleiotropy Test

#### 1. Simulate Pleiotropic Data with `pleiosim`

```r
library(pleiosim)

# Simulate pleiotropic genotypes and phenotypes
sim_data = pleiosim::pleiosim()
```

#### 2. Apply Pleiotropy Testing Methods

```r
# Run MixFisher with Davies approximation
result_mixfisher = run_mixfisher_davies(sim_data)
head(result_mixfisher)

# Run MixVar with Liu-modified approximation
result_mixvar = run_mixvar_liumod(sim_data)
head(result_mixvar)

# Run MinP test
result_minp = run_minp(sim_data)
head(result_minp)
```

#### 3. Run Multiple Pleiotropy Tests Using `pleiotest`

```r
# Run a comprehensive pleiotropy analysis
results = pleiotest(sim_data)

# View the first few results
head(results@pleiotest_stats)

# Check execution times of each test
results@execution_time
```

#### 4. Select Specific Tests

If you want to run only selected tests, you can specify them as follows:

```r
# Run only MixFisher and MinP tests
results_selected = pleiotest(sim_data, 
                             run_mixfisher_davies_test = TRUE, 
                             run_minp_test = TRUE, 
                             run_mixvar_liu_test = FALSE)

# View the output
head(results_selected@pleiotest_stats)
```

## Available Tests

`pleiotest` provides a range of statistical tests for detecting pleiotropy, including:



| Function                | Method               | Description                                                  |
| ----------------------- | -------------------- | ------------------------------------------------------------ |
| `run_amatz`             | AMATZ                | Adaptive multi-trait association test                        |
| `run_cmats`             | CMATS                | Cauchy-based multi-trait association test                    |
| `run_dsum`              | DSUM                 | Summation-based multi-trait pleiotropy test                  |
| `run_emats`             | EMATS                | Expectation-maximization multi-trait association test        |
| `run_metacca`           | metaCCA              | Canonical correlation analysis for meta-GWAS                 |
| `run_minp`              | MinP                 | Minimum p-value approach for pleiotropy detection            |
| `run_mixada`            | MixAda               | Adaptive weighted association method                         |
| `run_mixfisher_davies`  | MixFisher (Davies)   | Fisher’s method using Davies approximation                   |
| `run_mixfisher_liu`     | MixFisher (Liu)      | Fisher’s method using Liu approximation                      |
| `run_mixfisher_liumod`  | MixFisher (Liu-mod)  | Fisher’s method using Liu-modified approximation             |
| `run_mixsd_davies`      | MixSD (Davies)       | Standard deviation-based method using Davies approximation   |
| `run_mixsd_liu`         | MixSD (Liu)          | Standard deviation-based method using Liu approximation      |
| `run_mixsd_liumod`      | MixSD (Liu-mod)      | Standard deviation-based method using Liu-modified approximation |
| `run_mixtippett_davies` | MixTippett (Davies)  | Tippett’s combination test using Davies approximation        |
| `run_mixtippett_liu`    | MixTippett (Liu)     | Tippett’s combination test using Liu approximation           |
| `run_mixtippett_liumod` | MixTippett (Liu-mod) | Tippett’s combination test using Liu-modified approximation  |
| `run_mixvar_davies`     | MixVar (Davies)      | Variance-based pleiotropy test using Davies approximation    |
| `run_mixvar_liu`        | MixVar (Liu)         | Variance-based pleiotropy test using Liu approximation       |
| `run_mixvar_liumod`     | MixVar (Liu-mod)     | Variance-based pleiotropy test using Liu-modified approximation |
| `run_pcaq`              | PCAQ                 | Principal component adaptive quantile approach               |
| `run_pcfisher`          | PCFisher             | Principal component-based Fisher combination test            |
| `run_pclc`              | PCLC                 | Principal Component Likelihood Combination test              |
| `run_pcminp`            | PCMinP               | Minimum p-value across principal components                  |
| `run_pco`               | PCO                  | Principal Component Optimization method                      |
| `run_sum`               | SUM                  | Summation-based multi-trait pleiotropy test                  |
| `run_tates`             | TATES                | Trait-based association test that accounts for trait correlation |
| `run_vc`                | VC                   | Variance Component test for pleiotropy                       |
| `run_wald`              | Wald                 | Wald test for pleiotropy                                     |
| `run_wi`                | WI                   | Weighted integration of association signals                  |

## License

This package is licensed under the **[MIT License](LICENSE)**.

## Contributing

Contributions and suggestions are welcome! Feel free to open an [issue](https://github.com/your-github-username/pleiotest/issues) or submit a pull request.

## Contact

For questions or collaborations, contact **Wanjun Gu** at [wanjun.gu@ucsf.edu](mailto:wanjun.gu@ucsf.edu).


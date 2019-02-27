-   [Introduction](#introduction)
-   [Problem setting](#problem-setting)
-   [References](#references)

Introduction
============

This is a document to briefly explain the well-posedness condition in
Shapiro, Xie, and Zhang (2018) and provide an example of checking the
well-posedness condition.

Problem setting
===============

Consider the following Minimum Rank Matrix Completion (MRMC) problem:
$$
\\min\_{Y\\in \\mathbb R^{n\_1\\times n\_2}} {\\rm rank}(Y) {\\rm\\,\\,\\, subject\\, to\\,} Y\_{ij} =  M\_{ij}
$$
 **Definition** (Well-posedness condition).  
*We say that a matrix $Y\\in\\mathcal M\_r$ is well-posed for MRMC
problem, if *P*<sub>*Ω*</sub>(*Y*)=*M* and the following condition
holds*
$$
\\label{wp}
\\tag{1}
\\mathbb V\_{\\Omega^c}\\cap \\mathcal T\_{\\mathcal M\_r}(Y) = \\{0\\}.
$$

In Shapiro, Xie, and Zhang (2018), we also propose an equivalent
condition of well-posedness condition, which we can apply in practice.

**Theorem III.3**(Equivalent condition of well-posedness)

*Matrix $Y\\in \\mathcal M\_r$ satisfies condition ($\\ref{wp}$) if and
only if for any left side complement *F* and right side complement *G*
of *Y*, the column vectors
*g*<sub>*j*</sub><sup>⊤</sup> ⊗ *f*<sub>*i*</sub>,
(*i*, *j*)∈*Ω*<sup>*c*</sup>, are linearly independent.*

It is easy to verify that if there exists one left side complement *F*
and right side complement *G* the above condition is satisfied, then for
any left side complement *F* and right side complement *G* is satisfied.
According to it, the following is a function (in R) for checking
well-posedness condition.

    require(Matrix)
    require(rTensor)
    require(MASS)
    check_wp<-function(Pomega, Y){
      n1 = dim(Y)[1]
      n2 = dim(Y)[2]
      m = sum(!is.na(Pomega))
      VV = t(Null(Y))
      WW = Null(t(Y))
      WW_VV = kronecker_list(list(t(WW), VV))
      VEC_X = as.vector(Pomega)
      OmegaC = which(is.na(VEC_X))
      return(rankMatrix(WW_VV[,OmegaC])== n1*n2-m)
    }

**Example** 6 × 6 well-posed matrix (Wilson and Worcester 1939).

$$
M = \\left(
\\begin{array}{cccccc}
\*&0.56&0.16&0.48&0.24&0.64\\\\
0.56&\*&0.20&0.66 & 0.51&0.86\\\\
0.16&0.20&\*&0.18&0.07&0.23\\\\
0.48&0.66&0.18&\*&0.30&0.72\\\\
0.24&0.51&0.07&0.30&\*&0.41\\\\
0.64&0.86&0.23&0.72&0.41&\*
\\end{array}\\right)
$$
 In Wilson and Worcester (1939), it is proved there are only two
solutions of it, *D*1 = (0.64, 0.85, 0.06, 0.56, 0.50, 0.93) and
*D*2 = (0.42, 0.90, 0.06, 0.55, 0.39, 1.00). Therefore, the MRMC problem
should be well-posed since both solution are locally unique.

Now we check the well-posedness condition of these two solution.

    M = matrix(c(0,0.56,0.16,0.48,0.24,0.64,
                 0.56,0,0.2,0.66,0.51,0.86,
                 0.16,0.2, 0,0.18,0.07,0.23,
                 0.48,0.66,0.18,0,0.30,0.72,
                 0.24,0.51,0.07,0.3,0,0.41,
                 0.64,0.86,0.23,0.72,0.41,0
                 ), ncol = 6)
    Pomega = matrix(1, ncol = 6, nrow =6) - diag(NA,6,6)
    Y1 = M + diag(c(0.64,0.85,0.06,0.56,0.5,0.93))
    Y2 = M + diag(c(432/1015,1173/1300,311/4900,711/1300,116/300,998/1000))

Now, we can check whether both solutions are well-posed.

    check_wp(Pomega, Y1)

    ## [1] TRUE

    check_wp(Pomega, Y2)

    ## [1] TRUE

As we can see, to check well-posedness condition, we need to know the
true matrix. Therefore in Shapiro, Xie, and Zhang (2018), with a
numerical experiment, we show that when
$r \\leq \\mathfrak R(n\_1,n\_2,n\_1n\_2p) - 3$, $Y\\in\\mathcal M\_r$
and each entry is being observed with probability *p*, the MRMC problem
is well-posed with high probability. For a square matrix, we have the
following relationship between $\\frac{r+3}{n}$ and *p*:
$$
p\\geq \\frac{r+3}{n}(2-\\frac{r+3}{n})
$$

    x = seq(0,1,0.05)
    y = x*(2-x)
    plot(x,y, xlab = '(r+3)/n', ylab = 'p', type = 'l')

<img src="figure/rmd-plot-1.png" style="display: block; margin: auto;" />

------------------------------------------------------------------------

References
==========

Shapiro, Alexander, Yao Xie, and Rui Zhang. 2018. “Matrix Completion
with Deterministic Pattern - a Geometric Perspective,” February.
<https://arxiv.org/pdf/1802.00047>.

Wilson, Edwin B, and Jane Worcester. 1939. “The Resolution of Six Tests
into Three General Factors.” *Proceedings of the National Academy of
Sciences of the United States of America* 25 (2). National Academy of
Sciences: 73.

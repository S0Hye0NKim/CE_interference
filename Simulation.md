Causal Effect with interference
================
Sohyeon Kim

``` r
library(tidyverse)
```

# 1\. Simulation

This is a scenario for binary outcomes, heterogenoeus group between
groups.

## 1-1. Step 1

Of the 4000 individuals in the population,
![x\_{ij}](https://latex.codecogs.com/png.latex?x_%7Bij%7D "x_{ij}") is
generated as follows :

  
![x\_{ij} = \\begin{cases}0 &\\text{for randomly selected 480
samples}\\\\1 &\\text{for another 480 random samples}\\\\2&\\text{for
remaining individuals}
\\end{cases}](https://latex.codecogs.com/png.latex?x_%7Bij%7D%20%3D%20%5Cbegin%7Bcases%7D0%20%26%5Ctext%7Bfor%20randomly%20selected%20480%20samples%7D%5C%5C1%20%26%5Ctext%7Bfor%20another%20480%20random%20samples%7D%5C%5C2%26%5Ctext%7Bfor%20remaining%20individuals%7D%20%5Cend%7Bcases%7D
"x_{ij} = \\begin{cases}0 &\\text{for randomly selected 480 samples}\\\\1 &\\text{for another 480 random samples}\\\\2&\\text{for remaining individuals} \\end{cases}")  

``` r
set.seed(0)
x0 <- rep(0, 480)
x1 <- rep(1, 480)
x2 <- rep(2, 4000-480*2)

x <- sample(c(x0, x1, x2))
```

Potential outcome is defined as follows :

  
![y\_{ij}(z\_{ij}, \\alpha\_{ij})=\\begin{cases} g\_iz\_{ij}&\\text{if
}x\_{ij}=2\\text{ and }i=1,2\\\\g\_i(1-z\_{ij})&\\text{if
}x\_{ij}=2\\text{ and
}i=3,4\\\\x\_{ij}&\\text{otherwise}\\end{cases}](https://latex.codecogs.com/png.latex?y_%7Bij%7D%28z_%7Bij%7D%2C%20%5Calpha_%7Bij%7D%29%3D%5Cbegin%7Bcases%7D%20g_iz_%7Bij%7D%26%5Ctext%7Bif%20%7Dx_%7Bij%7D%3D2%5Ctext%7B%20and%20%7Di%3D1%2C2%5C%5Cg_i%281-z_%7Bij%7D%29%26%5Ctext%7Bif%20%7Dx_%7Bij%7D%3D2%5Ctext%7B%20and%20%7Di%3D3%2C4%5C%5Cx_%7Bij%7D%26%5Ctext%7Botherwise%7D%5Cend%7Bcases%7D
"y_{ij}(z_{ij}, \\alpha_{ij})=\\begin{cases} g_iz_{ij}&\\text{if }x_{ij}=2\\text{ and }i=1,2\\\\g_i(1-z_{ij})&\\text{if }x_{ij}=2\\text{ and }i=3,4\\\\x_{ij}&\\text{otherwise}\\end{cases}")  

## 1-2. Step 2

Groups were assigned ![\\alpha\_0,
\\alpha\_1](https://latex.codecogs.com/png.latex?%5Calpha_0%2C%20%5Calpha_1
"\\alpha_0, \\alpha_1") and individuals assigned z = 1,0 using two stage
permutation randomization.

``` r
g <- sample(c(1, 1, 0, 0))

z <- list()
for(i in 1:length(g)) {
  if(g[i] == 1) {
    z[[i]] <- sample(c(rep(0, 500), rep(1, 500)))
  } else {z[[i]] <- sample(c(rep(0, 800), rep(1, 200)))}
}

data <-  bind_cols(z) %>%
  `colnames<-`(1:4) %>%
  gather(key = "group_idx", value = z) %>%
  bind_cols(data.frame(g = rep(g, each = 1000), x = x))

y <- rep(NA, 4000)
for(i in 1:nrow(data)) {
  if(data$x[i] == 2) {
    y[i] <- ifelse(data$group_idx[i] %in% 1:2, data$g[i]*data$z[i], data$g[i]*(1-data$z[i]))
  } else {y[i] <- x[i]}
}

data <- mutate(data, y = y)
```

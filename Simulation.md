Causal Effect with interference
================
Sohyeon Kim

``` r
library(tidyverse)
library(pander)
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

## 1-3. Step 3

Calculate various causal effect estimators defined in the paper.

``` r
summary <- data %>%
  select(-x) %>%
  group_by(group_idx) %>%
  split(.$group_idx) %>%
  lapply(FUN = function(x) data.frame(Total_treated = sum(x$z), Cases_treated = sum(x$z*(x$y)), 
                                      Total_control = sum(1 - x$z), Cases_control = sum((1-x$z)*x$y))) %>%
  bind_rows(.id = "group_idx") %>%
  mutate(g_i = g) %>%
  select(group_idx, g_i, Total_treated, Cases_treated, Total_control, Cases_control)

summary %>% pander()
```

| group\_idx | g\_i | Total\_treated | Cases\_treated | Total\_control | Cases\_control |
| :--------: | :--: | :------------: | :------------: | :------------: | :------------: |
|     1      |  1   |      500       |      434       |      500       |       63       |
|     2      |  1   |      500       |      445       |      500       |       67       |
|     3      |  0   |      200       |       22       |      800       |       96       |
|     4      |  0   |      200       |       24       |      800       |       89       |

![\\text{Total treated }= \\sum\_jZ\_{ij} \\text{ when }Z\_{ij}
= 1](https://latex.codecogs.com/png.latex?%5Ctext%7BTotal%20treated%20%7D%3D%20%5Csum_jZ_%7Bij%7D%20%5Ctext%7B%20when%20%7DZ_%7Bij%7D%20%3D%201
"\\text{Total treated }= \\sum_jZ_{ij} \\text{ when }Z_{ij} = 1")

![\\text{Cases treated }= \\sum\_jZ\_{ij}Y\_{ij}(Z\_{i}) \\text{ when
}Z\_{ij}
= 1](https://latex.codecogs.com/png.latex?%5Ctext%7BCases%20treated%20%7D%3D%20%5Csum_jZ_%7Bij%7DY_%7Bij%7D%28Z_%7Bi%7D%29%20%5Ctext%7B%20when%20%7DZ_%7Bij%7D%20%3D%201
"\\text{Cases treated }= \\sum_jZ_{ij}Y_{ij}(Z_{i}) \\text{ when }Z_{ij} = 1")

![\\text{Total treated }= \\sum\_j(1-Z\_{ij}) \\text{ when }Z\_{ij}
= 1](https://latex.codecogs.com/png.latex?%5Ctext%7BTotal%20treated%20%7D%3D%20%5Csum_j%281-Z_%7Bij%7D%29%20%5Ctext%7B%20when%20%7DZ_%7Bij%7D%20%3D%201
"\\text{Total treated }= \\sum_j(1-Z_{ij}) \\text{ when }Z_{ij} = 1")

![\\text{Cases treated }= \\sum\_j(1-Z\_{ij})Y\_{ij}(Z\_{i}) \\text{
when }Z\_{ij}
= 1](https://latex.codecogs.com/png.latex?%5Ctext%7BCases%20treated%20%7D%3D%20%5Csum_j%281-Z_%7Bij%7D%29Y_%7Bij%7D%28Z_%7Bi%7D%29%20%5Ctext%7B%20when%20%7DZ_%7Bij%7D%20%3D%201
"\\text{Cases treated }= \\sum_j(1-Z_{ij})Y_{ij}(Z_{i}) \\text{ when }Z_{ij} = 1")

![g\_i=1](https://latex.codecogs.com/png.latex?g_i%3D1 "g_i=1") when
group assignment of group i is
![\\psi](https://latex.codecogs.com/png.latex?%5Cpsi "\\psi")

![g\_i=0](https://latex.codecogs.com/png.latex?g_i%3D0 "g_i=0") when
group assignment of group i is
![\\phi](https://latex.codecogs.com/png.latex?%5Cphi "\\phi")

### Direct Effect

``` r
summary %>% mutate(PO_1 = Cases_treated/Total_treated, 
                   PO_0 = Cases_control/Total_control) %>%
  select(group_idx, PO_1, PO_0) %>% 
  split(.$group_idx) %>%
  lapply(FUN = function(x) data.frame(CE_D = x$PO_0 - x$PO_1)) %>%
  bind_rows(.id = "grou_idx") %>%
  mutate(g = g) %>%
  group_by(g) %>%
  summarise(CE_D = mean(CE_D)) %>% pander
```

| g |  CE\_D   |
| :-: | :------: |
| 0 | 0.000625 |
| 1 | \-0.749  |

  
![\\bar{CE}^D(\\psi)= \\sum\_{i=1}^2\\frac{\\bar{CE}\_i^D(\\psi)}{2} =
-0.749](https://latex.codecogs.com/png.latex?%5Cbar%7BCE%7D%5ED%28%5Cpsi%29%3D%20%5Csum_%7Bi%3D1%7D%5E2%5Cfrac%7B%5Cbar%7BCE%7D_i%5ED%28%5Cpsi%29%7D%7B2%7D%20%3D%20-0.749
"\\bar{CE}^D(\\psi)= \\sum_{i=1}^2\\frac{\\bar{CE}_i^D(\\psi)}{2} = -0.749")  

  
![\\bar{CE}^D(\\phi)= \\sum\_{i=1}^2\\frac{\\bar{CE}\_i^D(\\phi)}{2}
= 0.00625](https://latex.codecogs.com/png.latex?%5Cbar%7BCE%7D%5ED%28%5Cphi%29%3D%20%5Csum_%7Bi%3D1%7D%5E2%5Cfrac%7B%5Cbar%7BCE%7D_i%5ED%28%5Cphi%29%7D%7B2%7D%20%3D%200.00625
"\\bar{CE}^D(\\phi)= \\sum_{i=1}^2\\frac{\\bar{CE}_i^D(\\phi)}{2} = 0.00625")  

### Indirect Effect

``` r
summary %>% mutate(PO_0 = Cases_control/Total_control) %>%
  group_by(g_i) %>%
  summarise(PO_0 = mean(PO_0)) %>% pander
```

| g\_i | PO\_0  |
| :--: | :----: |
|  0   | 0.1156 |
|  1   |  0.13  |

  
![\\begin{aligned}\\bar{CE}^I(\\phi,
\\psi)&=\\bar{Y}(0;\\phi)-\\bar{Y}(0;\\psi)\\\\&#10;
&=\\sum\_{i=3}^4\\frac{\\bar{Y}\_i(0;\\phi)}{2}-\\sum\_{i=1}^2\\frac{\\bar{Y}\_i(0;\\psi)}{2}\\\\&=0.116-0.13\\\\&=-0.0144\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cbar%7BCE%7D%5EI%28%5Cphi%2C%20%5Cpsi%29%26%3D%5Cbar%7BY%7D%280%3B%5Cphi%29-%5Cbar%7BY%7D%280%3B%5Cpsi%29%5C%5C%0A%20%20%20%20%26%3D%5Csum_%7Bi%3D3%7D%5E4%5Cfrac%7B%5Cbar%7BY%7D_i%280%3B%5Cphi%29%7D%7B2%7D-%5Csum_%7Bi%3D1%7D%5E2%5Cfrac%7B%5Cbar%7BY%7D_i%280%3B%5Cpsi%29%7D%7B2%7D%5C%5C%26%3D0.116-0.13%5C%5C%26%3D-0.0144%5Cend%7Baligned%7D
"\\begin{aligned}\\bar{CE}^I(\\phi, \\psi)&=\\bar{Y}(0;\\phi)-\\bar{Y}(0;\\psi)\\\\
    &=\\sum_{i=3}^4\\frac{\\bar{Y}_i(0;\\phi)}{2}-\\sum_{i=1}^2\\frac{\\bar{Y}_i(0;\\psi)}{2}\\\\&=0.116-0.13\\\\&=-0.0144\\end{aligned}")  

### Total Effect

``` r
summary %>% mutate(PO_1 = Cases_treated/Total_treated, 
                   PO_0 = Cases_control/Total_control) %>%
  group_by(g_i) %>%
  summarise(PO_1 = mean(PO_1), PO_0 = mean(PO_0)) %>% pander
```

| g\_i | PO\_1 | PO\_0  |
| :--: | :---: | :----: |
|  0   | 0.115 | 0.1156 |
|  1   | 0.879 |  0.13  |

  
![\\begin{aligned}\\bar{CE}^T(\\phi,\\psi)&=\\bar{Y}(0;\\phi)-\\bar{Y}(1;\\psi)\\\\&=\\sum\_{i=3}^4\\frac{\\bar{Y}\_i(0;\\phi)}{3}-\\sum\_{i=1}^2\\frac{\\bar{Y}\_i(1;\\psi)}{2}\\\\&=0.116-0.879=-0.763\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cbar%7BCE%7D%5ET%28%5Cphi%2C%5Cpsi%29%26%3D%5Cbar%7BY%7D%280%3B%5Cphi%29-%5Cbar%7BY%7D%281%3B%5Cpsi%29%5C%5C%26%3D%5Csum_%7Bi%3D3%7D%5E4%5Cfrac%7B%5Cbar%7BY%7D_i%280%3B%5Cphi%29%7D%7B3%7D-%5Csum_%7Bi%3D1%7D%5E2%5Cfrac%7B%5Cbar%7BY%7D_i%281%3B%5Cpsi%29%7D%7B2%7D%5C%5C%26%3D0.116-0.879%3D-0.763%5Cend%7Baligned%7D
"\\begin{aligned}\\bar{CE}^T(\\phi,\\psi)&=\\bar{Y}(0;\\phi)-\\bar{Y}(1;\\psi)\\\\&=\\sum_{i=3}^4\\frac{\\bar{Y}_i(0;\\phi)}{3}-\\sum_{i=1}^2\\frac{\\bar{Y}_i(1;\\psi)}{2}\\\\&=0.116-0.879=-0.763\\end{aligned}")  

  
![\\begin{aligned}\\bar{CE}^T(\\phi,\\psi)&=\\bar{CE}^I(\\phi,\\psi)+\\bar{CE}^D(\\psi)\\\\&=-0.0144-0.749\\\\&=-0.763\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cbar%7BCE%7D%5ET%28%5Cphi%2C%5Cpsi%29%26%3D%5Cbar%7BCE%7D%5EI%28%5Cphi%2C%5Cpsi%29%2B%5Cbar%7BCE%7D%5ED%28%5Cpsi%29%5C%5C%26%3D-0.0144-0.749%5C%5C%26%3D-0.763%5Cend%7Baligned%7D
"\\begin{aligned}\\bar{CE}^T(\\phi,\\psi)&=\\bar{CE}^I(\\phi,\\psi)+\\bar{CE}^D(\\psi)\\\\&=-0.0144-0.749\\\\&=-0.763\\end{aligned}")  

### Overall Effect

``` r
summary %>%
  mutate(PO = (Cases_treated + Cases_control)/(Total_treated + Total_control)) %>%
  group_by(g_i) %>%
  summarise(PO = mean(PO)) %>% pander
```

| g\_i |   PO   |
| :--: | :----: |
|  0   | 0.1155 |
|  1   | 0.5045 |

  
![\\begin{aligned}\\bar{CE}^O(\\phi,\\psi)&=\\bar{Y}(\\phi)-\\bar{Y}(\\psi)\\\\&=\\sum\_{i=3}^4\\frac{\\bar{Y}\_i(\\phi)}{2}-\\sum\_{i=1}^2\\frac{\\bar{Y}\_i(\\psi)}{2}\\\\&=0.116-0.505=-0.389\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cbar%7BCE%7D%5EO%28%5Cphi%2C%5Cpsi%29%26%3D%5Cbar%7BY%7D%28%5Cphi%29-%5Cbar%7BY%7D%28%5Cpsi%29%5C%5C%26%3D%5Csum_%7Bi%3D3%7D%5E4%5Cfrac%7B%5Cbar%7BY%7D_i%28%5Cphi%29%7D%7B2%7D-%5Csum_%7Bi%3D1%7D%5E2%5Cfrac%7B%5Cbar%7BY%7D_i%28%5Cpsi%29%7D%7B2%7D%5C%5C%26%3D0.116-0.505%3D-0.389%5Cend%7Baligned%7D
"\\begin{aligned}\\bar{CE}^O(\\phi,\\psi)&=\\bar{Y}(\\phi)-\\bar{Y}(\\psi)\\\\&=\\sum_{i=3}^4\\frac{\\bar{Y}_i(\\phi)}{2}-\\sum_{i=1}^2\\frac{\\bar{Y}_i(\\psi)}{2}\\\\&=0.116-0.505=-0.389\\end{aligned}")  

# Reference

Hudgens, M. G., & Halloran, M. E. (2008). Toward causal inference with
interference. Journal of the American Statistical Association, 103(482),
832-842.

Liu, L., & Hudgens, M. G. (2014). Large sample randomization inference
of causal effects in the presence of interference. Journal of the
American Statistical Association, 109(505), 288-301.




> Written with [StackEdit](https://stackedit.io/).

# Principal Component Analysis

## High dimentional data

 - The given Genedata is a high-dimentional data in which the , i.e, the given data contain multiple features. By reducing the given data to the lower dimentions it reveals the hidden simplified structuring of the data.

 PCA explores and calculate the maximum variance of the linear combinantions of the variables in the data and and follows the same for the second linear combinations and so on which results in the uncorrelation  of the factors.

**PCA** reduces the the dimantions of the given data by plotting/projecting them into the lower dimentiions geometrically which are also called the *Principal Components*

## Principal components
The first Principal Component reduces the distance between the Data and the PC while maximising the variance of the projections
The second principal component is chosen similarly to the firts *PC* and they are uncorrelated with the first PC

## PCA With Eigen Decomposition

Transforming the given data points X into the orthonormal matix P in **Z  = PX**
Where **X** contsins all the independent features of the data
```R
X = rbind( matrix(genedata$`x1`,1,301) , matrix(genedata$`x2`,1,301) , matrix(genedata$`x3`,1,301) ,
           matrix(genedata$`x4`,1,301) , matrix(genedata$`x4`,1,301) )
 ```



This can be done by removing the mean values in the given data matrix :

```R
mean_X_rows = matrix( rowMeans( X[1:5,] ) , 5 , 301 ) # Mean matrix of rowise means

X_features_de_mean = X[1:5,] - mean_X_rows # remove mean of the features from X. Note: features are from rows 1 to 5

```

The construction of the covariance Matrix is given by C<sub>x</sub>=(1/n)XX<sup>T</sup>
Where n = 301 as there are 301 Observations in the given Gene Data

```R
Cx = (1/301) * ( X_features_de_mean %*% t(X_features_de_mean) )
```

We can caluculate the eigen vectors and values of the Co-variance Matrix by ``` eigen()``` function

```R
﻿﻿﻿﻿﻿﻿﻿Eig = eigen(Cx)
Eig_values = Eig$values
Eig_vectors = Eig$vectors
```

```R
>print(Eig_values)
[1] 0.533426402 0.218550283 0.018060507 0.005359376 0.000000000
> print(Eig_vectors)
            [,1]        [,2]       [,3]        [,4]          [,5]
[1,] -0.37610477  0.42135554 -0.1162996  0.81699395  0.000000e+00
[2,] -0.09803589  0.84860104  0.2697961 -0.44438192  0.000000e+00
[3,] -0.73575378 -0.31968781  0.5902548 -0.08980714 -5.232037e-17
[4,] -0.39217378 -0.00823063 -0.5316407 -0.25197249 -7.071068e-01
[5,] -0.39217378 -0.00823063 -0.5316407 -0.25197249  7.071068e-01
```

Now if we want to transfoem the dimentions to K=2 dimentions then we select the eigen vectors of the covariance Matrix ***C<sub>x</sub>***  which are sorted in the decending values of the form **P**:

```R
k=2
P = Eig_vectors[1:k,]
```

This can be done by storing the eigen vectors in the P


Y=PX:(K x D)(D xN)=(K x N) :

```R
Y_k2 = P %*% X_features_de_mean # Y for k = 2
```

By ploting the reduced 2-Dimentional space Geometrically :
```R
plot( Y_k2[2,] , Y_k2[1,] , main="1,2") # for k = 2
```

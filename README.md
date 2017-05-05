# Linear-regression-sigma-70-promoters
Here multivariate linear regression has been used to predict the strength of sigma 70 core promoters.

Data set is 19 anderson promoters developed by anderson lab in UC Berkley - http://parts.igem.org/Promoters/Catalog/Anderson

The -10 and -35 hexamer motifs have been scored to predict the strength of sigma 70 promoters.
Their pssm was calculated and normalized as the two inputs.

Multivariant.py is the final code.

The optimised R-squared value(co-efficient of determination) was calculated to be 0.44 at a learning rate of 0.015 when gradient descent ran 1000000 times.

Another model will be generated soon  which doesn't use pssm, but frequence scores.

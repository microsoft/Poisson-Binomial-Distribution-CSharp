
# Poisson Binomial Distribution for C#

## About

Binomial probability distribution is the probability distribution of the sum of independent Bernoulli random variables with
non-uniform success probabilities. This implementation of the Poisson Binomial probability distribution is based on the python version developed by Mika Straka (https://github.com/tsakim/poibin).

Methods included:
* `pmf`: probability mass function
* `cdf`: cumulative distribution function

More information on the Poisson Binomial distribution can be found here: [Yili Hong, On computing the distribution function for the Poisson binomial distribution, Computational Statistics & Data Analysis, Volume 59, March 2013, pages 41-51,ISSN 0167-9473](http://dx.doi.org/10.1016/j.csda.2012.10.006)

## Usage

This code does not contain an FFT implementation. For it to run and work correctly, you will need to add a C# FFT library. I recommend and have used https://www.nayuki.io/res/free-small-fft-in-multiple-languages/FftTest.cs. If you use this source, the functional called in the code will call this FFT function without changes needed to the input parameters.

## Authors

Kay Toma, Carlos Zamora Cura

## Contributing

This project welcomes contributions and suggestions.  Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit https://cla.opensource.microsoft.com.

When you submit a pull request, a CLA bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

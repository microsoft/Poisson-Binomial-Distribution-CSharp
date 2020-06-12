    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using System.Numerics;
    using System.Diagnostics;
    using System.Globalization;
    using System.Collections;

    public class PoiBinProcessor : Processor
    {
    
        private bool debug = false;
        private bool debugNumerical = false;

        private double[] get_pmf_xi(double[] probabilities)
        {

            System.Numerics.Complex[] chi = new System.Numerics.Complex[probabilities.Length + 1];
            chi[0] = 1;
            int halfNumberTrials = (int)(probabilities.Length / 2 + probabilities.Length % 2);
            double omega = 2 * Math.PI / (probabilities.Length + 1);
            
            Complex[] firsthalfchi = get_chi(omega, probabilities, 1, halfNumberTrials + 1);

            // get the reverse order
            Complex[] secondHalfChi = new Complex[probabilities.Length - halfNumberTrials];

            int forward_index = 0;
            for (int p = secondHalfChi.Length - 1; p >= 0; p--)
            {
                secondHalfChi[forward_index] = firsthalfchi[p];
                forward_index++;
            }
                  
            // conjugate the second half of chis
            Complex[] conjugatedSecondHalfChi = new Complex[secondHalfChi.Length];
            for (int p = 0; p < secondHalfChi.Length; p++)
            {
                conjugatedSecondHalfChi[p] = Complex.Conjugate(secondHalfChi[p]);
            }

            // put the chi together and divide all chi by number_trials + 1
            int index;
            chi[0] = Complex.Divide(chi[0], new Complex(chi.Length, 0));
            for (index = 0; index < firsthalfchi.Length; index++)
            {
                if(debugNumerical)
                    Console.WriteLine("chi pre: " + firsthalfchi[index].ToString());
                chi[index + 1] = Complex.Divide(firsthalfchi[index], new Complex(chi.Length, 0));
                if(debugNumerical)
                    Console.WriteLine("chi post: " + chi[index + 1].ToString());
            }
            index++;
            for (int p = 0; p < secondHalfChi.Length; p++)
            {
                if(debugNumerical)
                    Console.WriteLine("chi pre: " + conjugatedSecondHalfChi[p].ToString());
                chi[index] = Complex.Divide(conjugatedSecondHalfChi[p], new Complex(chi.Length, 0));
                if(debugNumerical)
                    Console.WriteLine("chi post: " + chi[index].ToString());
                index++;
            }

            // FILL THIS IN WITH AN OPEN SOURCE FFT LIBRARY AND COMPUTE THE DISCRETE FOURIER TRANSFORM OF CHI (A COMPLEX VECTOR)
            Fft.Transform(chi, false);

            if (debugNumerical)
            {
                for (int f = 0; f < chi.Length; f++)
                {
                    Console.WriteLine("chi post fft: " + chi[f].ToString());
                }
            }

            double[] reals = new double[chi.Length];
            for (int p = 0; p < chi.Length; p++)
            {
                double real = chi[p].Real;
                reals[p] = real;
            }

            return reals;
            
        }

        private Complex[] get_chi(double omega, double[] success_probabilities, int start, int stop)
        {
            // create the idx_array
            // create the exp value per index of the array
            // create the new axis
            // exp_value = np.exp(self.omega * idx_array * 1j)
            // https://stackoverflow.com/questions/29241056/how-does-numpy-newaxis-work-and-when-to-use-it
            // https://docs.microsoft.com/en-us/dotnet/api/system.numerics.complex.exp?view=netframework-4.7.2
            if(debugNumerical)
                Console.WriteLine("Omega: " + omega.ToString());
            Complex i = new Complex(0, 1);
            int index = start;
            Complex[][] exp_value = new Complex[stop - start][];
            for (int p = 0; p < stop - start; p++)
            {
                Complex[] value = new Complex[] { Complex.Exp(omega * index * i) };
                exp_value[p] = value;
                if(debugNumerical)
                    Console.WriteLine("exp_value: " + value[0].ToString());
                index++;
            }
            
            // remainder of success probabilities
            // xy = 1 - self.success_probabilities + \
            Complex[] remainders = new Complex[success_probabilities.Length];
            for (int p = 0; p < success_probabilities.Length; p++)
            {
                remainders[p] = new Complex(1 - success_probabilities[p], 0);
            }

            // convert success_probabilities to complex type to be multiplied in next step
            Complex[] success_probabilities_complex = new Complex[success_probabilities.Length];
            for (int p = 0; p < success_probabilities.Length; p++)
            {
                success_probabilities_complex[p] = new Complex(success_probabilities[p], 0);
            }

            // self.success_probabilities * exp_value[:, np.newaxis]
            Complex[][] complexproducts = new Complex[exp_value.Length][];
            // cross multiple arrays
            // result should be an array of arrays
            // success probabilities: [1, 2, 3... n]
            // complex val: [[1], [2]...[n]]
            for (int p = 0; p < exp_value.Length; p++)
            {
                Complex[] complexvalue = exp_value[p];
                complexproducts[p] = multiply_arrays(success_probabilities_complex, complexvalue);
                // result of cross multiple of arrays of same length is to multiply 
                // [1, 1, 1] * [1, 2, 3] = [1, 2, 3]
            }

            // add the remainder to the product of the lists
            // success probabilities is a list
            // complex products is a list of lists where the interior lists are length of success_probabilities
            // xy = 1 - self.success_probabilities + \
            // self.success_probabilities * exp_value[:, np.newaxis]

            // have to be of the dimensions remainders = (x, ), complexproducts = (y, x)
            // add the single dimension array
            Complex[][] xy = new Complex[complexproducts.Length][];
            for (int p = 0; p < complexproducts.Length; p++)
            {

                if (complexproducts[p].Length != remainders.Length)
                {
                    // throw error
                    throw new System.ArgumentException("Length of complexproducts[p] does not equal length of remainders");
                }
                else
                {
                    Complex[] sum_list = new Complex[remainders.Length];
                    // remainders and complexproducts[p] should be same length;
                    for (int q = 0; q < complexproducts[p].Length; q++)
                    {
                        Complex sum = remainders[q] + complexproducts[p][q];
                        sum_list[q] = sum;
                    }
                    xy[p] = sum_list;
                }
            }

            // calculate arctan list of lists
            double[][] arctanresults = new double[xy.Length][];
            for (int p = 0; p < xy.Length; p++)
            {
                double[] arctanlist = new double[xy[p].Length];
                for (int q = 0; q < xy[p].Length; q++)
                {
                    double realcomponent = xy[p][q].Real;
                    double imaginarycomponent = xy[p][q].Imaginary;
                    double arctanresult = Math.Atan2(imaginarycomponent, realcomponent); // double check this is the right order and if arctan2 is same as numpy
                    arctanlist[q] = arctanresult;
                }
                arctanresults[p] = arctanlist;
            }


            // sum per list inside the list
            double[] argz_sum = new double[arctanresults.Length];
            for (int p = 0; p < arctanresults.Length; p++)
            {
                double sum = (double)0;
                for (int q = 0; q < arctanresults[p].Length; q++)
                {
                    sum += arctanresults[p][q];
                }
                argz_sum[p] = sum;
            }

            if (debugNumerical)
            {
                for (int k = 0; k < argz_sum.Length; k++)
                {
                    Console.WriteLine("argz_sum: " + argz_sum[k].ToString());
                }
            }
            // sum(axis=1):
            // array([[1, 1, 1],
            // [2, 2, 2],
            // [3, 3, 3]])
            // result --> [3, 6, 9]

            // calculate exponent arg
            // go through the list of lists and take the absolute value Complex.abs
            // natural log per element in list of lists
            // sum each of the lists
            double[] exparg = new double[xy.Length];
            for (int p = 0; p < xy.Length; p++)
            {
                double sum = (double)0;
                for (int q = 0; q < xy[p].Length; q++)
                {
                    Complex value = xy[p][q];
                    double abs = Complex.Abs(value);
                    double log = Math.Log(abs);
                    sum += log;
                }
                exparg[p] = sum;
            }
            // returns a list

            if (debugNumerical)
            {
                for (int e = 0; e < exparg.Length; e++)
                {
                    Console.WriteLine("exparg: " + exparg[e].ToString());
                }
            }

            // calculate d_value
            // for each value of exparg e^x
            double[] d_value = new double[exparg.Length];
            for (int p = 0; p < exparg.Length; p++)
            {
                double expArgValue = exparg[p];
                double dval = Math.Exp(expArgValue);

                d_value[p] = dval;
            }

            if (debugNumerical)
            {
                for (int d = 0; d < d_value.Length; d++)
                {
                    Console.WriteLine("d_value: " + d_value[d].ToString());
                }
            }

            // calculate chi
            // argz_sum * 1j
            // take the exponential of that product
            // multiple list by d_value
            Complex[] chi = new Complex[argz_sum.Length];
            for (int p = 0; p < argz_sum.Length; p++)
            {
                Complex argz_val = new Complex(argz_sum[p], 0);
                Complex product = argz_val * i;
                Complex expresult = Complex.Exp(product);

                double dval = d_value[p];
                Complex complexd_value = new Complex(dval, 0);
                chi[p] = complexd_value * expresult;
            }

            if (debugNumerical)
            {
                for (int c = 0; c < chi.Length; c++)
                {
                    Console.WriteLine("final chi: " + chi[c].ToString());
                }
            }

            return chi;
        }

        private Complex[] multiply_arrays(Complex[] success_probabilities, Complex[] value)
        {
            Complex[] product = new Complex[success_probabilities.Length];
            Complex complexvalue = value[0];
            for (int p = 0; p < success_probabilities.Length; p++)
            {
                Complex success = success_probabilities[p];
                product[p] = success * complexvalue;
            }

            return product;
        }

        private double[] get_cdf(double[] event_probabilities)
        {
            double[] cdf = new double[event_probabilities.Length];
            cdf[0] = event_probabilities[0];
            for (int p = 1; p < cdf.Length; p++)
            {
                cdf[p] = cdf[p - 1] + event_probabilities[p];
            }
            return cdf;
        }

    }

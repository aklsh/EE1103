# Random numbers

The aim of this assignment is to learn how to

* Use arguments (inputs) to main()
* Generate random numbers
* Estimate the Hamming distance between arrays of random bits,  without introducing any bit errors

Your program will have 3 inputs:

    N: the total number of bits
    M: an offset of the first location
    P: the number of random separation between bits, P < (N-M)

Create an array of N random bits of 0 and 1 (since N is an input to main, you will need to allocate/free memory).
Generate a random offset M, and store the bit value.
Generate P random numbers that are give you the separation of bits from M, and their bit values.
Use the Mth and the M+Pth bit values,  and compute the Hamming distance within the original array N by varying the offset (starting from the 0th location and going till the maximum possible location).
Check that you obtained the minimum Hamming distance for when the offset was M.

Check out this link for a clearer idea on hamming distance between two strings: https://www.geeksforgeeks.org/hamming-distance-two-strings/
.

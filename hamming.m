function w = hamming(n)
%HAMMING HAMMING(N) returns the N-point Hamming window.
w = .54 - .46*cos(2*pi*(0:n-1)'/(n-1));
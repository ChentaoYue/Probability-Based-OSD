function [OUT] = Combination(N,K)
%UNTITLED2 Summary of this function goes here
OUT = prod((N-K+1):N)/prod(1:K);
end


function [S SA] = Stirling2nd(n, k)
% S = Stirling2nd(N,K)
% N and K are integers
% Compute the Stirling's number of the second kind.
% It is the number of all possible partitions of the set {1:N}, where each
% (partition) has exactly K non-empty subsets.
%
% [S SA] = Stirling2nd(N,K) return a (N x K) array of all Stirling's
% numbers of the second kind of S(i,j) with i<=N, j<=min(K,i).
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History
%   Original: 17-May-2009
%   Last update: 18-May-2009, cosmetic changes
k=double(k);
n=double(n);
if k==0
    if n==0
        S = 1;
    else
        S = 0;
    end
    SA = zeros(n,0);
    return
end
SA = nan(n,k);
SA(:,1) = 1;
SA(sub2ind(size(SA),1:k,1:k)) = 1;
for i=2:n
    %for j=max(2,(i+k)-n):min(i-1,k) % ... % recursive path
    for j=2:min(i-1,k)
        SA(i,j) = SA(i-1,j-1) + j*SA(i-1,j);
    end
end
S = SA(n,k);
end
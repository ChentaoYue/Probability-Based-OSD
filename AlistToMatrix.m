function G = AlistToMatrix(N,K,D)
%% Transform the Alist format to Matrix
%% Chentao Yue. June, 2018
%% The university of Sydney, NSW

NonzeroColumn = ['AlistC_' num2str(N) '_' num2str(K) '_' num2str(D) '.txt'];
% NonzeroRow = ['AlistR_' num2str(N) '_' num2str(K) '_' num2str(D) '.txt'];

Columndata = load(NonzeroColumn);

H = zeros(N-K,N);
for i = 1:N
    Array = zeros(1,N-K);
    index = Columndata(i,:);
    index(index == 0) = [];
    Array(index) = 1;    
    H(:,i) = Array;   
end

[H2, P2] = GE(H);
H = H(:,P2);
P = H2(:,[(N-K+1):N]);
G = [ eye(K), P'];
H = [H(:,(N-K+1):N) H(:,1:(N-K))];

end

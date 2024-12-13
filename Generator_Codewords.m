function [c,G] = Generator_Codewords(p ,M)

N = p.N;
K = p.K;
D = p.D;
G = p.G;

c = mod(M*G,2);

end
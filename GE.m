function [G_1,P2,rank] = GE(G_1)

[K,N] = size(G_1);
P2 = 1:N;

i = 1;
j = 1;
offset = 1;
while (i <= K) && (j <= N)
    k = 0;
    tempcol = G_1(:,j);
    for ind = i:K
        if tempcol(ind)==1
            k = ind;
            break;
        end      
    end

   if k == 0
        check = K + offset;
        if check > N
           break; 
        end
        G_1(:,[j check]) = G_1(:,[check j]);
        P2([j check]) = P2([check j]);
        offset = offset+1;
        continue;  
   else   
        temp = G_1(k,:); 
        G_1 = abs(G_1 - tempcol*temp);   
        G_1(k,:) = G_1(i,:);
        G_1(i,:) = temp;
        i = i + 1;
        j = j + 1;
   end  
end
rank = i-1;

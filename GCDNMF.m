  function [X,G,cha_l,Q,H] = GCDNMF(A,C,k,alpha,beta)
% A is the adjacency matrix, C is the node attribute matrix
n = size(A,1); % number of nodes
m = size(C,2); % dimension of attributes
   
% random initialization    
G= rand(k,k);
Q = rand(n,k);   
U = rand(m,k);


%k-rank-D
H= K_rank_D_center(A,k);

if isnormalize
       H = H./max(realmin,repmat(sum(H,2),1,k)); %sum(H,2) 对H的每一行求和 返回一个列向量 repmat将列向量平铺为k个 组成一个n*k矩阵
end

cha = 1;
i = 0;
err=[];
cha_l=[];

while (cha > 1e-5&& i <100)  
    i = i+1;  
    
    % update G:
    G = G.*(H'*A*H + beta*Q'*Q)./max(realmin,(H'*H*G*H'*H + beta*G));   
    
    % update Q:
    Q = Q.*(alpha*C*U + 2*beta*Q*G)./max(realmin,(alpha*Q*U'*U + 2*beta*Q*Q'*Q));
    
    % update U:
   U = U.*(C'*Q)./max(realmin, (U*Q'*Q));
     
    % update H:
     H = H.*(2*A*H*G')./max(realmin,(2*H*G*H'*H*G'));
    
    % objective function Values
    first = norm(A-H*G*H','fro')^2; 
    second = alpha*norm(C-Q*U','fro')^2;
    third = beta*norm(G-Q'*Q,'fro')^2;    
    Ltemp = first + second + third;    
    
    err(size(err,1)+1,:) = Ltemp;

     if i==1
         cha=err(i);
     else
         cha=abs(err(i-1)-err(i));
     end
    cha_l(size(cha_l,1)+1,:)=cha;
    string=[num2str(i),' ', num2str(cha)];
    disp(string);

 end
%maximum assignment
   
 result1(1:n,1)=1:n;
 [~,index]=max(H,[],2);   %得到所属社团标号的列向量 将X中同一行的不同列的最大值提取出来 最大值结果存在~里，index里存的是每行最大值的列位置。
 result1(:,2)=index; 
 X=result1;
 

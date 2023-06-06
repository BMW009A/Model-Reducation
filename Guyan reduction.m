% Read mass matrix and stiffness matrix.
fileID_mass = fopen('your path');  %mass matrix
fileID_stiff = fopen('your path'); %stiffness matrix
Cmass = textscan(fileID_mass,'%f %f %f','Delimiter',',');
Cstiff = textscan(fileID_stiff,'%f %f %f','Delimiter',',');
N = 375;
fclose(fileID_mass);
fclose(fileID_stiff);
Cmtx = zeros(N, N);
Cstf = zeros(N, N);

for i = 1:length(Cmass{1,1})
    Cmtx(Cmass{1,1}(i),Cmass{1,2}(i)) = Cmass{1,3}(i);
end

for j = 1:length(Cstiff{1,1})
    Cstf(Cstiff{1,1}(j),Cstiff{1,2}(j)) = Cstiff{1,3}(j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Guyan Reduction Method %%%%%%%%%%%%%%%%%%%%%%%%%
wc = (max(sqrt(eig(Cstf,Cmtx))/(2*pi))*4) / (2 * pi)   % calculation max freq. of old mass martix and stiffness martix
max(sqrt(diag(Cstf)./diag(Cmtx)))

idx = [];
Index = [];
Value = [];

K = Cstf;
M = Cmtx;

while max(sqrt(diag(K)./diag(M))) > wc
    omg = sqrt(diag(K)./diag(M));
    % Put the index/value pair which value in diag bigger than freq into a table.
    for k = 1:length(omg)
        if omg(k) > wc
            Index = vertcat(Index,k);
            Value = vertcat(Value,omg(k));
        end
    end
    % Select the biggest value in the table but with the smallest index as
    % the partition of the matrix.
    MaxDig = table(Index, Value);
    [m_1, m_2] = max(MaxDig.Value);
    n = MaxDig.Index(m_2);
    if n == 1
        M11 = M(n+1:end,n+1:end);
        M12 = M(n+1:end,1:n);
        M21 = M(1:n,n+1:end);
        M22 = M(1:n,1:n);

        K11 = K(n+1:end,n+1:end);
        K12 = K(n+1:end,1:n);
        K21 = K(1:n,n+1:end);
        K22 = K(1:n,1:n);
     
    elseif n == length(omg)
        M11 = M(1:n-1,1:n-1);
        M12 = M(1:n-1,n:end);
        M21 = M(n:end,1:n-1);
        M22 = M(n:end,n:end);

        K11 = K(1:n-1,1:n-1);
        K12 = K(1:n-1,n:end);
        K21 = K(n:end,1:n-1);
        K22 = K(n:end,n:end);
        
    else
        M11 = [M(1:n-1,1:n-1) M(1:n-1,n+1:end); M(n+1:end,1:n-1) M(n+1:end,n+1:end)];
        M12 = [M(1:n-1,n); M(n+1:end,n)];
        M21 = [M(n,1:n-1) M(n,n+1:end)];
        M22 = [M(n,n)];

        K11 = [K(1:n-1,1:n-1) K(1:n-1,n+1:end); K(n+1:end,1:n-1) K(n+1:end,n+1:end)];
        K12 = [K(1:n-1,n); K(n+1:end,n)];
        K21 = [K(n,1:n-1) K(n,n+1:end)];
        K22 = [K(n,n)];
    end
    
    M = M11 - K12*pinv(K22)*M21-M12*pinv(K22)*K21+K12*pinv(K22)*M22*pinv(K22)*K21;
    K = K11 - K12*pinv(K22)*K21;

    Index = [];
    Value = [];
end

%%%%%%%%%%%% Error Calculation %%%%%%%%%%%%
% ORG = eig(Cstf,Cmtx);
% RED = eig(K,M);
% ORG_SORT = sort(ORG);
% RED_SORT = sort(RED);
% 
% err = (RED_SORT(1:75)- ORG_SORT(1:75))./ORG_SORT(1:75)
% [err_max, err_max_index] = max(err)
% [err_min, err_min_index] = min(err)

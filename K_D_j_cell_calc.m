function [ D_j_cell , K_j_cell , I_j_not]= K_D_j_cell_calc ( D_bar , K_bar, N , T , L)

D_j_cell= {};
K_j_cell= {};
I_j_not = {};
for j = 1 : N
    I_j = (j-1) * T +1 : j*T ;
    I_j_not{end+1} = setdiff( 1: N * T *(L+1) , I_j ) ; %*(L+1)
    D_j_cell{end+1} = D_bar ( I_j_not { end } , I_j_not { end });
    K_j_cell{end+1} = K_bar ( : , I_j_not { end });
end


end


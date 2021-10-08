function [K_bar , D_bar] = Kernel_ADDM_Network_topology ( Y , L , N , T , dataSize , param)
% K_bar = - ones ( L+1 , N , T , T ,1 ) ; % , dataSize
% Y size = T * N %* dataSize
for outer_counter = 1: L+1
    for inner_counter = 1: N
        K_bar_i_l = Kernel_i_l_ADDM_Network_topology ( squeeze ( Y (:, inner_counter) ), ...
            inner_counter , outer_counter-1 , T , dataSize , param);
        if (outer_counter == 1 && inner_counter==1)
            K_bar = K_bar_i_l; 
            D_bar = K_bar_i_l;
        else
            K_bar = horzcat ( K_bar ,  K_bar_i_l );
            D_bar = blkdiag ( D_bar , K_bar_i_l );
        end
              
    end
%     inner_D_bar = bdiag( squeeze ( K_bar( outer_counter , :, :,:, :) ) );
end

end


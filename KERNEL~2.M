function K_i_l = Kernel_i_l_ADDM_Network_topology ( y , i, l, T , dataSize, param )

% Y : T * 1
% y = squeeze ( Y(:, 1));
[A,B] = meshgrid(y,y);
K_i_l = [A(:) B(:)];
K_i_l = reshape ( K_i_l , T , T , 2);
% K_i_l = exp ( - (K_i_l(:,:,1)-K_i_l(:,:,2)).^2 / (2*param ));
K_i_l = K_i_l(:,:,1) .* K_i_l(:,:,2);
K_i_l = abs ( K_i_l ) + K_i_l.^2;


end


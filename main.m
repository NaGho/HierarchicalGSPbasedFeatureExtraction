%% In the Name of God


close all;
% clc();
%% GSP toolbox initialization
addpath(genpath ( 'C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\research\Graph Signal Processing\Brain Connectivity Matlab Toolbox'));
addpath(genpath ( 'C:\Users\Nafise\Dropbox\PhD\research\Graph Signal Processing\Brain Connectivity Matlab Toolbox'));
% cd ('C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\research\Graph Signal Processing\gspbox-0.7.4\gspbox')
% % cd ('C:\Users\Nafise\Dropbox\PhD\research\Graph Signal Processing\gspbox-0.7.4\gspbox')
% gsp_start
% gsp_make
% gsp_install
% % gsp_install_unlocbox
% % cd 3rdparty/sources/flann-1.8.4-src/cmake/
% % cmake .. 
% % make
% cd ('C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\research\Graph Signal Processing\matlab code')
% % cd ('C:\Users\Nafise\Dropbox\PhD\research\Graph Signal Processing\matlab code')
%% data loading
data_counter_preictal_array = 100 * ones ( 5,1);%[24 , 42 , 72, 97,30];  
for start_patient_counter = 1:5
    version_save = start_patient_counter + 10;
    start_data_counter = 1 ;%1
    patient_dog_name = 'Dog'; %'Patient'
    data_counter_preictal = data_counter_preictal_array(start_patient_counter);
    data_counter_interictal = data_counter_preictal ;
    patient_counter = 1;
    where_data = 'C:\1My files\Kaggle Data Set Prediction\' ;
%     where_data = 'E:\academic files\PhD\kaggle data set\' ;
    desired_time_horizon = 30; %3000000
    [interictal_DataSize_Graph_top , preictal_DataSize, whole_interictal_Data, train_interictal_time_length,train_interictal_sampling_frequency,channel_names, ...
        train_interictal_DataSize_Graph_top, valid_interictal_DataSize_Graph_top,test_interictal_DataSize, ...
       whole_preictal_Data, train_preictal_DataSize,  valid_preictal_DataSize, test_preictal_DataSize  ] = ...
        Data_Loader_Wrapped ('Dog' ,data_counter_interictal, data_counter_preictal, patient_counter,...
        start_data_counter , start_patient_counter , desired_time_horizon , where_data ); 
    [~ , ~ , time_horizon] = size ( whole_interictal_Data ) ; %5 ; 
    channel_num_total = length ( channel_names); % 4 ; 
    channel_locations = channel_locations_calc ( channel_num_total , patient_dog_name) ; 
    %% ADDM for Network topology Identification
    N = channel_num_total;
    T = time_horizon;
    Y_interictal = permute ( whole_interictal_Data , [3, 2 , 1] );% T , N , trainDataSize
    Y_interictal = Y_interictal (1:T,1:N ,:);
    fprintf ( 'ADDM graph learning Interictal Feature Extraction ... \n');
    ADDM_feature_Names = ['mean_Weights','Network_Density', 'clustering_coef_bd' ,  'efficiency_bin' , 'num_connected_components' ,  'size_largest_component',  ...
         ' average_num_neighbors', 'num_self_loops', 'Characteristic path length', 'global efficiency' , 'eccentricity', ...
         'radius', 'diameter', 'graph_weights'];
     L = 2; % memory in time
     % (N,N,L+1)
    [ ADDM_weights_temp_space_interictal , ADDM_mean_Weights_feature ,ADDM_Network_Density_feature ,ADDM_clustering_coef_bd_feature, ADDM_efficiency_bin_features, ...
        num_connected_components_feature , size_largest_component_feature , average_num_neighbors_feature ,num_self_loops_feature, ...
            lambda_feature , efficiency_feature , ecc_feature , radius_feature , diameter_feature ,  ADDM_All_Norm_Weights_features_interictal ] ...
            = ADDM_Network_topology_Identification ( N , T , Y_interictal , interictal_DataSize_Graph_top, 'interictal' , L);
    interictal_features_Graph_Learning = horzcat ( ADDM_mean_Weights_feature ,ADDM_Network_Density_feature ,ADDM_clustering_coef_bd_feature, ADDM_efficiency_bin_features, ...
        num_connected_components_feature , size_largest_component_feature , average_num_neighbors_feature ,num_self_loops_feature, ...
            lambda_feature , efficiency_feature , ecc_feature , radius_feature , diameter_feature );

    folder_to_save =strcat( strcat( 'C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\Eclipse Workspace\Prediction_Project\ADDM_features_v' , ...
        num2str(version_save)) , '\ADDM_features_');
 
%     folder_to_save = ...
%         'C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\eclipse-workspace-parallerl\Code_Graph_Project\ADDM_features_v13\ADDM_features_';
    csvwrite(strcat( folder_to_save , 'interictal_train.csv'),interictal_features_Graph_Learning(1:train_interictal_DataSize_Graph_top,:))
    csvwrite(strcat( folder_to_save ,'interictal_valid.csv'),interictal_features_Graph_Learning...
        ( train_interictal_DataSize_Graph_top + 1 : train_interictal_DataSize_Graph_top + valid_interictal_DataSize_Graph_top , : ))
    csvwrite(strcat( folder_to_save ,'interictal_test.csv'),interictal_features_Graph_Learning...
        ( train_interictal_DataSize_Graph_top + valid_interictal_DataSize_Graph_top +1:interictal_DataSize_Graph_top,:))



    Y_preictal = permute ( whole_preictal_Data , [3, 2 , 1] );% T , N , trainDataSize
    Y_preictal = Y_preictal (1:T,1:N ,:);
    fprintf ( 'ADDM graph learning Preictal Feature Extraction ... \n');
    [ ADDM_weights_temp_space_preictal , ADDM_mean_Weights_feature ,ADDM_Network_Density_feature ,ADDM_clustering_coef_bd_feature, ADDM_efficiency_bin_features, ...
        num_connected_components_feature , size_largest_component_feature , average_num_neighbors_feature ,num_self_loops_feature, ...
            lambda_feature , efficiency_feature , ecc_feature , radius_feature , diameter_feature , ADDM_All_Norm_Weights_features_preictal ] ...
            = ADDM_Network_topology_Identification ( N , T , Y_preictal , preictal_DataSize, 'preictal' , L);
    preictal_features_Graph_Learning = horzcat ( ADDM_mean_Weights_feature ,ADDM_Network_Density_feature ,ADDM_clustering_coef_bd_feature, ADDM_efficiency_bin_features, ...
        num_connected_components_feature , size_largest_component_feature , average_num_neighbors_feature ,num_self_loops_feature, ...
            lambda_feature , efficiency_feature , ecc_feature , radius_feature , diameter_feature );

    csvwrite(strcat( folder_to_save , 'preictal_train.csv'),preictal_features_Graph_Learning(1:train_preictal_DataSize,:))
    csvwrite(strcat( folder_to_save ,'preictal_valid.csv'),preictal_features_Graph_Learning...
        ( train_preictal_DataSize + 1 : train_preictal_DataSize + valid_preictal_DataSize,:))
    csvwrite(strcat( folder_to_save ,'preictal_test.csv'),preictal_features_Graph_Learning...
        (train_preictal_DataSize + valid_preictal_DataSize +1:preictal_DataSize,:))

    %% Strong Graph construction Spatial using mean of learned graph 
    GSP_num_features = 4;
    where_data = 'C:\1My files\Kaggle Data Set Prediction\' ;
    desired_time_horizon = T; %3000000
    data_counter_preictal = data_counter_preictal_array(start_patient_counter);
    data_counter_interictal = data_counter_preictal ;
    [interictal_DataSize , preictal_DataSize, whole_interictal_Data, train_interictal_time_length,train_interictal_sampling_frequency,channel_names, ...
        train_interictal_DataSize, valid_interictal_DataSize,test_interictal_DataSize, ...
       whole_preictal_Data, train_preictal_DataSize,  valid_preictal_DataSize, test_preictal_DataSize  ] = ...
        Data_Loader_Wrapped ('Dog' ,data_counter_interictal, data_counter_preictal, patient_counter,...
        start_data_counter , start_patient_counter , desired_time_horizon , where_data ); 
    [~ , ~ , time_horizon] = size ( whole_interictal_Data ) ; %5 ; 
    
    Mean_weights_diag_temp_space_interictal = squeeze(mean( ADDM_weights_temp_space_interictal(:,:,:,1),1)); 
    W_total = zeros ( N*T , N*T );     
    for T_counter = 1:T
        start_index = (T_counter-1)*N + 1;
        end_index = T_counter*N ;
        W_total(start_index:end_index, start_index:end_index) = Mean_weights_diag_temp_space_interictal;
        for L_counter = 2:L+1
            if ( end_index + ( L_counter - 1 ) * N <= N*T )
                Mean_weights_temp_space_interictal = squeeze(mean( ADDM_weights_temp_space_interictal(:,:,:,L_counter),1)); 
                % Mean_weights_temp_space_preictal 
                W_total ( start_index  : end_index ,...
                    end_index + (L_counter-2)*N +1 : end_index + (L_counter-1) * N ) = transpose (Mean_weights_temp_space_interictal );

                W_total ( end_index + (L_counter-2)*N +1 : end_index + (L_counter-1) * N  , ...
                    start_index : end_index ) = Mean_weights_temp_space_interictal;
            end
        end
    end
    
    
    G_total.W = W_total - diag(diag(W_total));
    G_total = gsp_graph_default_parameters(G_total);
    G_total.coords = zeros ( N * T , 3);
    G_total = gsp_compute_fourier_basis(G_total);
    EigenVectors_reshaped = reshape ( G_total.U , N , T , N * T);
    
    Add_largest_indices = 3; % 1:only values or 2: only indices, 3:both values and indices 4:values and indices for wavelet
    
    Y_interictal = permute ( whole_interictal_Data , [2, 3 , 1] );%  N , T , DataSize
    Y_interictal = reshape ( Y_interictal , time_horizon * channel_num_total , interictal_DataSize );%?????? careful while reshaping??????
    fprintf ( 'GSP Interictal Feature Extraction ... \n');
    [Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ...
        spectral_wavelet_feature, variation_feature , GSP_feature_names ] = GSP_based_features ( G_total , Y_interictal , 'interictal' ,Add_largest_indices, GSP_num_features);
    interictal_features_GSP = horzcat (Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ... %, Normalized_Graph_Laplacian_feature , random_walk_feature
        spectral_wavelet_feature, variation_feature);

%     folder_to_save = 'C:\1My files\Kaggle Data Set Prediction\GSP features\GSP_features_ ';   
    folder_to_save = ...
        strcat(strcat('C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\Eclipse Workspace\Prediction_Project\GSP_features_v' ,...
            num2str(version_save)),'\GSP_features_');
    csvwrite(strcat( folder_to_save , 'interictal_train.csv'),interictal_features_GSP(1:train_interictal_DataSize,:,:))
    csvwrite(strcat( folder_to_save ,'interictal_valid.csv'),interictal_features_GSP...
        ( train_interictal_DataSize + 1 : train_interictal_DataSize + valid_interictal_DataSize , : , : ))
    csvwrite(strcat( folder_to_save ,'interictal_test.csv'),interictal_features_GSP...
        ( train_interictal_DataSize + valid_interictal_DataSize +1:interictal_DataSize,:,:))




    Y_preictal = permute ( whole_preictal_Data , [2, 3 , 1] );%  N , T, DataSize
    Y_preictal = reshape ( Y_preictal , time_horizon * channel_num_total , preictal_DataSize );%?????? careful while reshaping??????
    fprintf ( 'GSP Preictal Feature Extraction ... \n');
    [Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ...
        spectral_wavelet_feature, variation_feature, GSP_feature_names ] = GSP_based_features ( G_total , Y_preictal , 'preictal' , Add_largest_indices , GSP_num_features);
    preictal_features_GSP = horzcat (Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ... %
        spectral_wavelet_feature, variation_feature);

    csvwrite(strcat( folder_to_save , 'preictal_train.csv'),preictal_features_GSP(1:train_preictal_DataSize,:,:))
    csvwrite(strcat( folder_to_save ,'preictal_valid.csv'),preictal_features_GSP...
        ( train_preictal_DataSize + 1 : train_preictal_DataSize + valid_preictal_DataSize,:,:))
    csvwrite(strcat( folder_to_save ,'preictal_test.csv'),preictal_features_GSP...
        (train_preictal_DataSize + valid_preictal_DataSize +1:preictal_DataSize,:,:))
end
%% show scattering plot of features
close all;
which = 'GSP'; % ADDM , GSP
if(0)
    folder_scatter = strcat(strcat(strcat('C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\eclipse-workspace-parallerl\Code_Graph_Project\',which),...
        '_features_v4\'), which);
else
    folder_scatter = strcat(strcat(strcat('C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\eclipse-workspace-parallerl\Code_Graph_Project\features_no_indices\',which),...
        '_features_v4\'), which);
end
if ( strcmp(which , 'ADDM'))
    which_features_to_see = 1:13;
else
    which_features_to_see = [1,89,97];
end
preictal_feature_scatter_plot = vertcat( csvread(strcat(folder_scatter, '_features_preictal_train.csv')) , ...
    csvread(strcat(folder_scatter, '_features_preictal_valid.csv')) , ...
    csvread(strcat(folder_scatter, '_features_preictal_test.csv') ));
interictal_feature_scatter_plot = vertcat( csvread( strcat(folder_scatter, '_features_interictal_train.csv')) , ...
    csvread(strcat(folder_scatter, '_features_interictal_valid.csv')) ,...
    csvread(strcat(folder_scatter, '_features_interictal_test.csv')) );
preictal_feature_scatter_plot = preictal_feature_scatter_plot (:,which_features_to_see);
interictal_feature_scatter_plot = interictal_feature_scatter_plot (:,which_features_to_see);
for counter1 = 1: length(which_features_to_see)
    for counter2 = counter1 + 1: length(which_features_to_see)
        figure();
        scatter( preictal_feature_scatter_plot(:,counter1), preictal_feature_scatter_plot(:,counter2))
        hold on;
        scatter( interictal_feature_scatter_plot(:,counter1) , interictal_feature_scatter_plot(:,counter2))
    end
end



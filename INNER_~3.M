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

% 0:mean of all channels and all time, 1: network density
outer_layer_input_option = 0 ; 

% 0: seperate inner-outer, all weightes leaned, 1: merged inner-outer, all weightes leaned
GSP_construction_option = 1 ; 

for start_patient_counter = [1,3,4,5]
    version_save = start_patient_counter + 10 ;
    start_data_counter = 1 ;%1
    patient_dog_name = 'Dog'; %'Patient'
    data_counter_preictal = data_counter_preictal_array(start_patient_counter);
    data_counter_interictal = data_counter_preictal ;
    patient_counter = 1;
    where_data = 'C:\1My files\Kaggle Data Set Prediction\' ;
%     where_data = 'E:\academic files\PhD\kaggle data set\' ;
    desired_time_horizon = 70; %3000000
    
    [interictal_DataSize_Graph_top , preictal_DataSize_Graph_top, whole_interictal_Data, train_interictal_time_length,train_interictal_sampling_frequency,channel_names, ...
        train_interictal_DataSize_Graph_top, valid_interictal_DataSize_Graph_top,test_interictal_DataSize_Graph_top, ...
       whole_preictal_Data, train_preictal_DataSize_Graph_top,  valid_preictal_DataSize_Graph_top, test_preictal_DataSize_Graph_top  ] = ...
        Data_Loader_Wrapped ('Dog' ,data_counter_interictal, data_counter_preictal, patient_counter,...
        start_data_counter , start_patient_counter , desired_time_horizon , where_data ); 
    [~ , ~ , time_horizon] = size ( whole_interictal_Data ) ; %5 ; 
    channel_num_total = length ( channel_names); % 4 ; 
    channel_locations = channel_locations_calc ( channel_num_total , patient_dog_name) ; 

    %% ADDM for Network topology Identification
    outer_temporal_resolution = 5;
    
    fprintf ( 'ADDM graph learning Interictal Feature Extraction ... \n');
    ADDM_feature_Names = ['mean_Weights','Network_Density', 'clustering_coef_bd' ,  'efficiency_bin' , 'num_connected_components' ,  'size_largest_component',  ...
         ' average_num_neighbors', 'num_self_loops', 'Characteristic path length', 'global efficiency' , 'eccentricity', ...
         'radius', 'diameter', 'graph_weights'];
     L = 1; % memory in time
     N = channel_num_total;
     preictal_features_Graph_Learning = [];
     interictal_features_Graph_Learning = [];
     outer_temporal_interictal_data = -1000 * ones (outer_temporal_resolution,1,interictal_DataSize_Graph_top);
     outer_temporal_preictal_data = -1000 * ones (outer_temporal_resolution,1,preictal_DataSize_Graph_top);
     ADDM_weights_temp_space_interictal_inner = -1000 * ones (outer_temporal_resolution,interictal_DataSize_Graph_top,N,N,L+1,T);
     ADDM_weights_temp_space_preictal_inner = -1000 * ones (outer_temporal_resolution,preictal_DataSize_Graph_top,N,N,L+1,T);
     for outer_temporal_counter = 1: outer_temporal_resolution
         % (N,N,L+1)
        T =  floor(desired_time_horizon/outer_temporal_resolution);
        Y_interictal = permute ( whole_interictal_Data ( : , 1 : N ,( outer_temporal_counter - 1 ) * T + 1 : outer_temporal_counter * T) , [3, 2 , 1] );% T , N , trainDataSize
        Y_interictal = Y_interictal (1:T,1:N ,:);
        [ ADDM_weights_temp_space_interictal_inner_temp , ADDM_mean_Weights_feature_interictal ,ADDM_Network_Density_feature_interictal ,ADDM_clustering_coef_bd_feature_interictal, ADDM_efficiency_bin_features_interictal, ...
            num_connected_components_feature_interictal , size_largest_component_feature_interictal , ...
            average_num_neighbors_feature_interictal ,num_self_loops_feature_interictal, ...
                lambda_feature_interictal , efficiency_feature_interictal , ecc_feature_interictal , ...
                radius_feature_interictal , diameter_feature_interictal ,  ADDM_All_Norm_Weights_features_interictal ] ...
                = ADDM_Network_topology_Identification ( N , T , Y_interictal , interictal_DataSize_Graph_top, 'interictal' , L);
         ADDM_weights_temp_space_interictal_inner(outer_temporal_counter,:,:,:,:,:) = reshape( ADDM_weights_temp_space_interictal_inner_temp, ...
             size(ADDM_weights_temp_space_interictal_inner(outer_temporal_counter,:,:,:,:,:)));
        interictal_features_Graph_Learning = horzcat ( interictal_features_Graph_Learning ,...
            horzcat ( ADDM_mean_Weights_feature_interictal ,ADDM_Network_Density_feature_interictal ,...
            ADDM_clustering_coef_bd_feature_interictal, ADDM_efficiency_bin_features_interictal, ...
            num_connected_components_feature_interictal , size_largest_component_feature_interictal ,...
            average_num_neighbors_feature_interictal ,num_self_loops_feature_interictal, ...
                lambda_feature_interictal , efficiency_feature_interictal , ecc_feature_interictal ,...
                radius_feature_interictal , diameter_feature_interictal ));
        
        
        Y_preictal = permute ( whole_preictal_Data( : , 1 : N ,( outer_temporal_counter - 1 ) * T + 1 : outer_temporal_counter * T) , [3, 2 , 1] );% T , N , trainDataSize
        Y_preictal = Y_preictal (1:T,1:N ,:);
        fprintf ( 'ADDM graph learning Preictal Feature Extraction ... \n');
        [ ADDM_weights_temp_space_preictal_inner_temp  , ADDM_mean_Weights_feature_preictal ,...
            ADDM_Network_Density_feature_preictal ,ADDM_clustering_coef_bd_feature_preictal, ADDM_efficiency_bin_features_preictal, ...
            num_connected_components_feature_preictal , size_largest_component_feature_preictal , ...
            average_num_neighbors_feature_preictal ,num_self_loops_feature_preictal, ...
                lambda_feature_preictal , efficiency_feature_preictal , ecc_feature_preictal , ...
                radius_feature_preictal , diameter_feature_preictal , ADDM_All_Norm_Weights_features_preictal ] ...
                = ADDM_Network_topology_Identification ( N , T , Y_preictal , preictal_DataSize_Graph_top, 'preictal' , L);
         ADDM_weights_temp_space_preictal_inner(outer_temporal_counter,:,:,:,:,:) = reshape( ADDM_weights_temp_space_preictal_inner_temp, ...
             size(ADDM_weights_temp_space_preictal_inner(outer_temporal_counter,:,:,:,:,:)));
         
        preictal_features_Graph_Learning = horzcat ( preictal_features_Graph_Learning , ...
            horzcat ( ADDM_mean_Weights_feature_preictal ,ADDM_Network_Density_feature_preictal ,ADDM_clustering_coef_bd_feature_preictal, ADDM_efficiency_bin_features_preictal, ...
            num_connected_components_feature_preictal , size_largest_component_feature_preictal , ...
            average_num_neighbors_feature_preictal ,num_self_loops_feature_preictal, ...
                lambda_feature_preictal , efficiency_feature_preictal , ecc_feature_preictal , radius_feature_preictal , diameter_feature_preictal ) );
        
        if(outer_layer_input_option == 0)
            outer_temporal_interictal_data(outer_temporal_counter,:,:) = squeeze ( mean(mean(Y_interictal, 2),1) );
            outer_temporal_preictal_data(outer_temporal_counter,:,:) = squeeze ( mean(mean(Y_preictal, 2),1) );
        elseif(outer_layer_input_option == 1)
            outer_temporal_interictal_data(outer_temporal_counter,:,:) = ADDM_Network_Density_feature_interictal;
            outer_temporal_preictal_data(outer_temporal_counter,:,:) = ADDM_All_Norm_Weights_features_preictal;
        end
       
     end
    N = 1;
    T =  outer_temporal_resolution;
    L = outer_temporal_resolution - 1 ; % memory in time
    Y_interictal = outer_temporal_interictal_data;
    
    [ ADDM_weights_temp_space_interictal_outer , ADDM_mean_Weights_feature ,ADDM_Network_Density_feature ,ADDM_clustering_coef_bd_feature, ADDM_efficiency_bin_features, ...
        num_connected_components_feature , size_largest_component_feature , average_num_neighbors_feature ,num_self_loops_feature, ...
            lambda_feature , efficiency_feature , ecc_feature , radius_feature , diameter_feature ,  ADDM_All_Norm_Weights_features_interictal ] ...
            = ADDM_Network_topology_Identification ( N , T , Y_interictal , interictal_DataSize_Graph_top, 'interictal' , L);
    interictal_features_Graph_Learning = horzcat ( interictal_features_Graph_Learning , horzcat ( ADDM_mean_Weights_feature ,ADDM_Network_Density_feature ,ADDM_clustering_coef_bd_feature, ADDM_efficiency_bin_features, ...
        num_connected_components_feature , size_largest_component_feature , average_num_neighbors_feature ,num_self_loops_feature, ...
            lambda_feature , efficiency_feature , ecc_feature , radius_feature , diameter_feature ));


    Y_preictal = outer_temporal_preictal_data;
    fprintf ( 'ADDM graph learning Preictal Feature Extraction ... \n');
    [ ADDM_weights_temp_space_preictal_outer , ADDM_mean_Weights_feature ,ADDM_Network_Density_feature ,ADDM_clustering_coef_bd_feature, ADDM_efficiency_bin_features, ...
        num_connected_components_feature , size_largest_component_feature , average_num_neighbors_feature ,num_self_loops_feature, ...
            lambda_feature , efficiency_feature , ecc_feature , radius_feature , diameter_feature , ADDM_All_Norm_Weights_features_preictal ] ...
            = ADDM_Network_topology_Identification ( N , T , Y_preictal , preictal_DataSize_Graph_top, 'preictal' , L);
    preictal_features_Graph_Learning = horzcat ( preictal_features_Graph_Learning , horzcat ( ADDM_mean_Weights_feature ,ADDM_Network_Density_feature ,ADDM_clustering_coef_bd_feature, ADDM_efficiency_bin_features, ...
        num_connected_components_feature , size_largest_component_feature , average_num_neighbors_feature ,num_self_loops_feature, ...
            lambda_feature , efficiency_feature , ecc_feature , radius_feature , diameter_feature ) );
            
            
     folder_to_save =strcat( strcat( 'C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\Eclipse Workspace\Prediction_Project\ADDM_features_v' , ...
    num2str(version_save)) , '\ADDM_features_');
%     folder_to_save = ...
%         'C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\eclipse-workspace-parallerl\Code_Graph_Project\ADDM_features_v13\ADDM_features_';
    csvwrite(strcat( folder_to_save , 'interictal_train.csv'),interictal_features_Graph_Learning(1:train_interictal_DataSize_Graph_top,:))
    csvwrite(strcat( folder_to_save ,'interictal_valid.csv'),interictal_features_Graph_Learning...
        ( train_interictal_DataSize_Graph_top + 1 : train_interictal_DataSize_Graph_top + valid_interictal_DataSize_Graph_top , : ))
    csvwrite(strcat( folder_to_save ,'interictal_test.csv'),interictal_features_Graph_Learning...
        ( train_interictal_DataSize_Graph_top + valid_interictal_DataSize_Graph_top +1:interictal_DataSize_Graph_top,:))

    csvwrite(strcat( folder_to_save , 'preictal_train.csv'),preictal_features_Graph_Learning(1:train_preictal_DataSize_Graph_top,:))
    csvwrite(strcat( folder_to_save ,'preictal_valid.csv'),preictal_features_Graph_Learning...
        ( train_preictal_DataSize_Graph_top + 1 : train_preictal_DataSize_Graph_top + valid_preictal_DataSize_Graph_top,:))
    csvwrite(strcat( folder_to_save ,'preictal_test.csv'),preictal_features_Graph_Learning...
        (train_preictal_DataSize_Graph_top + valid_preictal_DataSize_Graph_top +1:preictal_DataSize_Graph_top,:))   
        
    %% Strong Graph construction Spatial using mean of learned graph 
    ADDM_weights_temp_space_interictal_inner = real(ADDM_weights_temp_space_interictal_inner);
    ADDM_weights_temp_space_interictal_outer = real(ADDM_weights_temp_space_interictal_outer);
    GSP_num_features = 4;
    if(GSP_construction_option == 2)
        desired_time_horizon = 200; %3000000   
        [interictal_DataSize_GSP , preictal_DataSize_GSP, whole_interictal_Data, train_interictal_time_length,train_interictal_sampling_frequency,channel_names, ...
            train_interictal_DataSize_GSP, valid_interictal_DataSize_GSP,test_interictal_DataSize_GSP, ...
           whole_preictal_Data, train_preictal_DataSize_GSP,  valid_preictal_DataSize_GSP, test_preictal_DataSize_GSP  ] = ...
            Data_Loader_Wrapped ('Dog' ,data_counter_interictal, data_counter_preictal, patient_counter,...
            start_data_counter , start_patient_counter , desired_time_horizon , where_data ); 
        [~ , ~ , time_horizon] = size ( whole_interictal_Data ) ; %5 ; 
        channel_num_total = length ( channel_names); % 4 ; 
        channel_locations = channel_locations_calc ( channel_num_total , patient_dog_name) ; 
        
    elseif(GSP_construction_option == 0 || GSP_construction_option == 1)
        
        interictal_DataSize_GSP = interictal_DataSize_Graph_top;
        preictal_DataSize_GSP = preictal_DataSize_Graph_top;
        train_interictal_DataSize_GSP = train_interictal_DataSize_Graph_top;
        valid_interictal_DataSize_GSP = valid_interictal_DataSize_Graph_top;
        test_interictal_DataSize_GSP = test_interictal_DataSize_Graph_top;
        train_preictal_DataSize_GSP = train_preictal_DataSize_Graph_top;
        valid_preictal_DataSize_GSP = valid_preictal_DataSize_Graph_top;
        test_preictal_DataSize_GSP = test_preictal_DataSize_Graph_top;    
    end   
    interictal_features_GSP = [];
    preictal_features_GSP = [];
    if(GSP_construction_option == 0 )
        for outer_temporal_counter = 1: outer_temporal_resolution
             % (N,N,L+1)
            N = channel_num_total;
            T =  floor(desired_time_horizon/outer_temporal_resolution);
            W_total = weight_Calculator(squeeze(ADDM_weights_temp_space_interictal_inner(outer_temporal_counter,:,:,:,:)),N,T,L);

            G_total.W = W_total - diag(diag(W_total));
            G_total = gsp_graph_default_parameters(G_total);
            G_total.coords = zeros ( N * T , 3);
            G_total = gsp_compute_fourier_basis(G_total);
            EigenVectors_reshaped = reshape ( G_total.U , N , T , N * T);

            Add_largest_indices = 3; % 1:only values or 2: only indices, 3:both values and indices 4:values and indices for wavelet


            fprintf ( 'GSP Interictal Feature Extraction ... \n');
            Y_interictal = permute ( whole_interictal_Data ( : , 1 : N ,( outer_temporal_counter - 1 ) * T + 1 : outer_temporal_counter * T) , [3, 2 , 1] );% T , N , trainDataSize
            Y_interictal = Y_interictal (1:T,1:N ,:);
            Y_interictal = reshape ( Y_interictal , T * N , size(Y_interictal,3) );

            [Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ...
                spectral_wavelet_feature, variation_feature , GSP_feature_names ] = ...
                GSP_based_features ( G_total , Y_interictal , 'interictal' ,Add_largest_indices, GSP_num_features );
            interictal_features_GSP = horzcat ( interictal_features_GSP , horzcat (Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ... %, Normalized_Graph_Laplacian_feature , random_walk_feature
                spectral_wavelet_feature, variation_feature));


            Y_preictal = permute ( whole_preictal_Data  ( : , 1 : N ,( outer_temporal_counter - 1 ) * T + 1 : outer_temporal_counter * T), [2, 3 , 1] );%  N , T, DataSize
            Y_preictal = reshape ( Y_preictal , T * N , size(Y_preictal,3) );%?????? careful while reshaping??????
            fprintf ( 'GSP Preictal Feature Extraction ... \n');
            [Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ...
                spectral_wavelet_feature, variation_feature, GSP_feature_names ] = ...
                GSP_based_features ( G_total , Y_preictal , 'preictal' , Add_largest_indices , GSP_num_features );
            preictal_features_GSP = horzcat ( preictal_features_GSP , horzcat (Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ... %
                spectral_wavelet_feature, variation_feature));   
        end

        N = 1;
        T =  outer_temporal_resolution;
        W_total = weight_Calculator(ADDM_weights_temp_space_interictal_outer,N,T,L);
        G_total.W = W_total - diag(diag(W_total));
        G_total = gsp_graph_default_parameters(G_total);
        G_total.coords = zeros ( N * T , 3);
        G_total = gsp_compute_fourier_basis(G_total);
        EigenVectors_reshaped = reshape ( G_total.U , N , T , N * T);    
        Add_largest_indices = 3; % 1:only values or 2: only indices, 3:both values and indices 4:values and indices for wavelet

        Y_interictal = outer_temporal_interictal_data;
        Y_interictal = reshape ( Y_interictal , T * N , size(Y_interictal,3) );
        fprintf ( 'GSP Interictal Feature Extraction ... \n');
        [Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ...
            spectral_wavelet_feature, variation_feature , GSP_feature_names ] = ...
            GSP_based_features ( G_total , Y_interictal , 'interictal' ,Add_largest_indices , GSP_num_features );
        interictal_features_GSP = horzcat ( interictal_features_GSP , horzcat (Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ... %, Normalized_Graph_Laplacian_feature , random_walk_feature
            spectral_wavelet_feature, variation_feature));

        Y_preictal = outer_temporal_preictal_data ; 
        Y_preictal = reshape ( Y_preictal , T * N , size(Y_preictal,3) );%?????? careful while reshaping??????
            fprintf ( 'GSP Preictal Feature Extraction ... \n');
        [Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ...
            spectral_wavelet_feature, variation_feature, GSP_feature_names ] = ...
            GSP_based_features ( G_total , Y_preictal , 'preictal' , Add_largest_indices , GSP_num_features );
        preictal_features_GSP = horzcat ( preictal_features_GSP , horzcat (Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ... %
            spectral_wavelet_feature, variation_feature) );
    end
    if ( GSP_construction_option == 1 )
        
        N = channel_num_total;
        T =  time_horizon;
        
        W_total = weight_Calculator_merged_inner_outer(ADDM_weights_temp_space_interictal_inner,ADDM_weights_temp_space_interictal_outer);
        G_total.W = W_total - diag(diag(W_total));
        G_total = gsp_graph_default_parameters(G_total);
        G_total.coords = zeros ( N * T , 3);
        G_total = gsp_compute_fourier_basis(G_total);
        EigenVectors_reshaped = reshape ( G_total.U , N , T , N * T);    
        Add_largest_indices = 3; % 1:only values or 2: only indices, 3:both values and indices 4:values and indices for wavelet

        Y_interictal = permute ( whole_interictal_Data , [3, 2 , 1] )  ;
        Y_interictal = reshape ( Y_interictal , T * N , size(Y_interictal,3) );
        fprintf ( 'GSP Interictal Feature Extraction ... \n');
        [Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ...
            spectral_wavelet_feature, variation_feature , GSP_feature_names ] = ...
            GSP_based_features ( G_total , Y_interictal , 'interictal' ,Add_largest_indices , GSP_num_features );
        interictal_features_GSP = horzcat ( interictal_features_GSP , horzcat (Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ... %, Normalized_Graph_Laplacian_feature , random_walk_feature
            spectral_wavelet_feature, variation_feature));

        Y_preictal = permute ( whole_preictal_Data , [3, 2 , 1] )  ; 
        Y_preictal = reshape ( Y_preictal , T * N , size(Y_preictal,3) );%?????? careful while reshaping??????
            fprintf ( 'GSP Preictal Feature Extraction ... \n');
        [Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ...
            spectral_wavelet_feature, variation_feature, GSP_feature_names ] = ...
            GSP_based_features ( G_total , Y_preictal , 'preictal' , Add_largest_indices , GSP_num_features );
        preictal_features_GSP = horzcat ( preictal_features_GSP , horzcat (Graph_Laplacian_feature, Normalized_Graph_Laplacian_feature , random_walk_feature, ... %
            spectral_wavelet_feature, variation_feature) );
        
        
    end
    
    
%     folder_to_save = 'C:\1My files\Kaggle Data Set Prediction\GSP features\GSP_features_ ';   
    folder_to_save = ...
        strcat(strcat('C:\Users\Nafiseh Ghoroghchian\Dropbox\PhD\Eclipse Workspace\Prediction_Project\GSP_features_v' ,...
            num2str(version_save)),'\GSP_features_');
        
        
    csvwrite(strcat( folder_to_save , 'interictal_train.csv'),interictal_features_Graph_Learning(1:train_interictal_DataSize_GSP,:))
    csvwrite(strcat( folder_to_save ,'interictal_valid.csv'),interictal_features_Graph_Learning...
        ( train_interictal_DataSize_GSP + 1 : train_interictal_DataSize_GSP + valid_interictal_DataSize_GSP , : ))
    csvwrite(strcat( folder_to_save ,'interictal_test.csv'),interictal_features_Graph_Learning...
        ( train_interictal_DataSize_GSP + valid_interictal_DataSize_GSP +1:interictal_DataSize_GSP,:))

    csvwrite(strcat( folder_to_save , 'preictal_train.csv'),preictal_features_Graph_Learning(1:train_preictal_DataSize_GSP,:))
    csvwrite(strcat( folder_to_save ,'preictal_valid.csv'),preictal_features_Graph_Learning...
        ( train_preictal_DataSize_GSP + 1 : train_preictal_DataSize_GSP + valid_preictal_DataSize_GSP,:))
    csvwrite(strcat( folder_to_save ,'preictal_test.csv'),preictal_features_Graph_Learning...
        (train_preictal_DataSize_GSP + valid_preictal_DataSize_GSP +1:preictal_DataSize_GSP,:))   
    
end

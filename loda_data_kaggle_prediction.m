function [ data, time_length,sampling_frequency,channels ]...
    = loda_data_kaggle_prediction( where_data , type_data , data_counter , patient_counter , patiante_dog_name...
    , sampling_precedure , sampling_parameter , start_data_counter ,  start_patient_counter)
% sampling_precedure can be 'normal_undersampling',  'window_mean' ,
% '1_difference'
% sampling_parameter for 'normal_undersampling' is the desired time-horizon
% sampling_parameter for 'window_mean' is the desired time-horizon
% 
address = get_address_loda_data_kaggle_prediction ( 1, start_patient_counter ,  type_data , where_data , patiante_dog_name);
% data a matrix of EEG sample values arranged row x column as electrode x time
general_struct = load(address);
general_struct = general_struct.(genvarname(strcat ( strcat(type_data , '_segment_') , num2str(1) )));
% data_length_sec: the time duration of each data row
time_length = general_struct.data_length_sec ;
% sampling_frequency: the number of data samples representing 1 second of EEG data
sampling_frequency = general_struct.sampling_frequency ;
% channels: a list of electrode names corresponding to the rows in the data field
channels = general_struct.channels;
% sequence: the index of the data segment within the one hour series of clips. For example, preictal_segment_6.mat 
% has a sequence number of 6, and represents the iEEG data from 50 to 60 minutes into the preictal data.

sampling_indices = ceil ( linspace (1 , time_length * sampling_frequency ,  sampling_parameter)); % 1:2:100


if ( strcmp (sampling_precedure , 'normal_undersampling') )
    counter_gen = 1;
    for outer_counter = start_patient_counter : start_patient_counter + patient_counter - 1
        for  inner_counter = start_data_counter:data_counter+start_data_counter - 1
            address = get_address_loda_data_kaggle_prediction ( inner_counter, outer_counter ,...
                type_data , where_data , patiante_dog_name );
            % data a matrix of EEG sample values arranged row x column as electrode x time
            if exist(address, 'file') == 2
                general_struct = load(address);
                general_struct = general_struct.(genvarname(strcat ( strcat(type_data , '_segment_') , num2str(inner_counter) )));
                data (counter_gen,:,:) = reshape ( general_struct.data(:,sampling_indices) ,...
                    1, length ( channels) ,length(sampling_indices) ) ; %undersamplinggg??????????
            else
                break;
            end
            counter_gen = counter_gen + 1;
            
        end
    end
    
    
    
elseif ( strcmp (sampling_precedure , 'window_mean') )
    window_duration = sampling_indices(2) - sampling_indices(1);
    counter_gen = 1;
    for outer_counter = start_patient_counter : start_patient_counter + patient_counter - 1
        for  inner_counter = start_data_counter : data_counter + start_data_counter - 1
            address = get_address_loda_data_kaggle_prediction ...
                ( inner_counter, outer_counter ,  type_data , where_data , patiante_dog_name);
            % data a matrix of EEG sample values arranged row x column as electrode x time
            if ( exist(address, 'file') == 2 )
                general_struct = load(address);
                general_struct = general_struct.(genvarname(strcat ( strcat(type_data , '_segment_') , num2str(inner_counter) )));
                moving_avg = movmean( general_struct.data, window_duration ,2 );
                if( size(moving_avg(:,sampling_indices),1) * size(moving_avg(:,sampling_indices),2) == ...
                    length ( channels) * length(sampling_indices) ) 
                    data (counter_gen,:,:) = reshape ( moving_avg(:,sampling_indices)  ,...
                        1, length ( channels) ,length(sampling_indices) ) ; %undersamplinggg??????????
                end
            else
                break;
            end
            counter_gen = counter_gen + 1;
%             fprintf('inner counter  = %d\n', inner_counter/data_counter )
        end
    end
    
    
elseif ( strcmp (sampling_precedure , '1_difference') )
    
    counter_gen = 1;
    for outer_counter = start_patient_counter : start_patient_counter + patient_counter - 1
        for  inner_counter = start_data_counter : data_counter + start_data_counter - 1
            address = get_address_loda_data_kaggle_prediction ...
                ( inner_counter, outer_counter ,  type_data , where_data , patiante_dog_name);
            % data a matrix of EEG sample values arranged row x column as electrode x time
            if exist(address, 'file') == 2
                general_struct = load(address);
                general_struct = general_struct.(genvarname(strcat ( strcat(type_data , '_segment_') , num2str(inner_counter) )));
                difference_1_signal = diff( general_struct.data , 2 );
                sign_diff_indices = find([0 diff(sign(difference_1_signal))]~=0);
                if (length(sign_diff_indices)>sampling_parameter)
                    sampling_indices = ceil ( linspace (1 , time_length * sampling_frequency ,  sampling_parameter));
                else
                    sampling_indices = sign_diff_indices;
                end
                data (counter_gen,:,:) = reshape ( difference_1_signal(:,sampling_indices)  ,...
                    1, length ( channels) ,length(sampling_indices) ) ; %undersamplinggg??????????
            else
                break;
            end
            counter_gen = counter_gen + 1;
        end
    end
    
      
end





end


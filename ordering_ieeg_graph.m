function F_total = ordering_ieeg_graph ( F_ieeg )

[channel_num_total , time_horizon] = size (F_ieeg);
F_total = reshape( F_ieeg,[ channel_num_total*time_horizon , 1 ]);


end


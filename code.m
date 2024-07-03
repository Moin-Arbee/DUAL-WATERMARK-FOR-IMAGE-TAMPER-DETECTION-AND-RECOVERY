% Load the MAT file
load('image.mat');

% Access the variable 'image'
input_image = image;
%disp(input_image);
%imshow(input_image);
M = 64;
%input image is of size 64 * 64


% Initialize variables
blocks = M/2; % Number of blocks in each dimension

block_size = [2, 2];
image_size = [64, 64];

N = M/2 * M/2 ; %total number of 2 x 2 blocks in the input_image.

% Initialize matrix to store X(i,j) values for each block
X = zeros(blocks, blocks);

% Calculate X(i,j) for each block
for i = 1:blocks
    for j = 1:blocks
        X(i,j) = ((i-1) * (M/2)) + (j-1);
    end
end


% Initialize prime number k
k = 17;

% Initialise the polynomial 
% Define the divisor polynomial in binary

divisor_binary = '111';

%the only irreducible degree-2 polynomial is x^2 +x+1.

% Initialize matrix to store X' values for each block
X_prime = zeros(blocks, blocks);


% Display the dimensions of X_prime matrix
%
% 
% disp(['Dimensions of X_prime matrix: ', num2str(size(X_prime))]);


% Calculate X'(i,j) for each block
for i = 1:blocks
    for j = 1:blocks
        X_prime(i,j) = mod(k * X(i,j), N) + 1;
    end
end



% Display X_prime matrix
% disp('X_prime:');
%disp(X_prime);





% 1. Fill the original number array from X-prime matrix, by traversing its first row
original_numbers = X_prime(1, :);
% Print original_numbers
%disp('Original numbers:');
%disp(original_numbers);

% 2. Get the permuted order
rng(k); % Set the random number generator seed
permuted_order = randperm(length(original_numbers));

% 3. Store it in a permutation map
permutation_map = containers.Map(original_numbers, original_numbers(permuted_order));

% Display the permuted numbers
% disp('Original -> Permuted:');
%disp([original_numbers; cell2mat(values(permutation_map))]);


% 4. Create a matrix, X-prime-prime, of size 32x32
X_prime_prime = zeros(32, 32);

% 5. Traverse X-prime's first row, find the corresponding mapping of it in permutation map
for i = 1:32
    original_number = X_prime(1, i);
    permuted_number = permutation_map(original_number);
    
    % 6. Find the position of the value obtained via map in X-prime
    % Find that column number, corresponding to the permuted number in X-prime
    [~, column_number] = find(X_prime == permuted_number, 1, 'first'); 

    % Fill the column of X-prime-prime with the column of original_number from X-prime
    X_prime_prime(:, column_number) = X_prime(:, i);
end





% mapping process completed
%The original image is first divided horizontally into two equal parts

% The original image is first divided horizontally into two equal parts 
% Blocks A and B are partner-block which are located at the same position 
% in these two parts. They are both of
% size 2×2 pixels.

% we are partitioning at after 16th rows.
% all entries in 17th row and first row follow this property - 
% that |M1 − M2| = N/2 our N = 1024.





watermarked_image = input_image;

watermark_matrix = cell(M/4, M/2);
%we will traverse the input image in blocks of 2 x 2, upper(a) and lower
%half( partner of a, i.e. b)
% then, we check the corresponding value in X-prime-prime for both blocks
% x and y


positions = {[0, 0], [0, 1], [1, 0], [1, 1]};

   for row = 1:2:31  %1:2:31
         for col = 1:2:63  %1:2:63
             % Initialize sum variable for the 2x2 block
            block_sum = 0;
            for pos = 1:length(positions)
                i = positions{pos}(1);
                j = positions{pos}(2);
                x = row + i;
                y = col + j;
                 % Add pixel value to block_sum
               
                curr_pixel_val = input_image(x,y);
                
                 % Extract the last three bits of the pixel value
                last_three_bits = bitand(curr_pixel_val, 7); % Equivalent to temp_pixel_value AND 7 in binary
                
                % Convert the last three bits to decimal value
                last_three_bits_decimal = bin2dec(dec2bin(last_three_bits));
                
                % Subtract the last three bits from the pixel value
                temp_pixel_value_without_last_three_bits = curr_pixel_val - last_three_bits_decimal;

                input_image(x,y)=temp_pixel_value_without_last_three_bits;
               
                block_sum = block_sum + input_image(x, y);
            end
            % Compute average pixel value for the block
            average_value_a = block_sum / 4;
            
            % to do the same for its partner in the lower half
            % it has same columns but it has it's row = row + 32
            partner_row = row + 32;
            block_sum = 0;
            for pos = 1:length(positions)
                i = positions{pos}(1);
                j = positions{pos}(2);
                x = partner_row + i;
                y = col + j;
                 % Add pixel value to block_sum
                curr_pixel_val = input_image(x,y);
                 % Extract the last three bits of the pixel value
                last_three_bits = bitand(curr_pixel_val, 7); % Equivalent to temp_pixel_value AND 7 in binary
                
                % Convert the last three bits to decimal value
                last_three_bits_decimal = bin2dec(dec2bin(last_three_bits));
                
                % Subtract the last three bits from the pixel value
                temp_pixel_value_without_last_three_bits = curr_pixel_val - last_three_bits_decimal;

                input_image(x,y)=temp_pixel_value_without_last_three_bits;
                block_sum = block_sum + input_image(x, y);
            end
            % Compute average pixel value for the block
            average_value_b= block_sum / 4;
            
            % now, we have average intensities of 2 x 2 partner blocks a
            % and b
            % now, to calculate p and v, USING CRC.
            binary_value_a = dec2bin(round(average_value_a), 8); % Convert to 8-bit binary
            binary_value_b = dec2bin(round(average_value_b), 8); % Convert to 8-bit binary
            

            % Extract starting 5 MSBs from both values
            msb_a = binary_value_a(1:5);
            msb_b = binary_value_b(1:5);
            
            % Convert MSBs back to strings
            msb_string_a = num2str(msb_a);
            msb_string_b = num2str(msb_b);
            
            % Concatenate MSB strings
            msb_concatenated = strcat(msb_string_a, msb_string_b);
            
            % Concatenate additional 2 bits for CRC
            dividend = strcat(msb_concatenated, '00');
            %disp("Dividend ");
           % disp(dividend);
            % divisor is defined at line 34.
            
            % Perform CRC using mod2div function
            %disp("Remainder ");
            remainder = mod2div(dividend, divisor_binary);
            %disp(remainder);
            % Extract the last two bits of the remainder
            to_append = remainder(end-1:end);

            
            % Replace the last two bits of the dividend with the remainder
            % disp("Hello");
            dividend(end-1:end) = to_append;
          
            
            % Calculate the watermark information
            watermark_info = dividend;
           % display(watermark_info);
            %disp(watermark_info);
            % Store the watermark information in the watermark_matrix
            watermark_matrix{(row+1)/2, (col+1)/2} = watermark_info;
            

        end
    end

%we have stored the 12bit information for original input image in
%watermark_matrix
%Now, we embed this information into original image, referring to the X''
  
% 2.1.3. Block watermark embedding

% smoothing input_image first
% Initialize the smoothed_input_image matrix

%disp(w_m = watermark_matrix{1, 2});
block_num = 1;
smoothed_input_image = input_image;

% Iterate over each value in the watermark_matrix (16 x 32)
 for row = 1:16
     for col = 1:32
        w_m = watermark_matrix{row, col};
        binary_string = w_m;

        [row_num, col_num] = blockNumberToRowCol(X_prime_prime(row, col), block_size, image_size);
     
        %we have to embed the info in row_num and col_num in the original
        %image
        start = 1;
        ending = 3;
        
        for pos = 1:4
                x = positions{pos}(1);
                y = positions{pos}(2);
                i = row_num + x;
                j = col_num + y;
               
                % Get the value of that pixel in smoothed_input_image 
                x_val = smoothed_input_image(i, j);
               
                current_bits = binary_string(start:ending);
                start = start + 3;
                ending = ending + 3;
              
                decimal_value = bin2dec(current_bits);
            
                % Extract the 3 least significant bits of X
                X_LSB = bitand(x_val, 7); % Equivalent to X AND (111) in binary
                
                % Calculate the difference between the watermark information and the LSBs
                v = decimal_value - X_LSB;
                %disp(v);
                
                % Apply the smoothing function
                if abs(v) < 5
                    smoothed_input_image(i, j) = x_val - X_LSB + v;
                elseif v >= 5
                    smoothed_input_image(i, j) = x_val - X_LSB + v + 8;
                else
                    smoothed_input_image(i, j) = x_val - X_LSB + v - 8;
                end
                
                % embed in smoothed_input_image
                 % Extract the current pixel value from smoothed_input_image
                temp_pixel_value = x_val;
                
                % Extract the last three bits of the pixel value
                last_three_bits = bitand(temp_pixel_value, 7); % Equivalent to temp_pixel_value AND 7 in binary
                
                % Convert the last three bits to decimal value
                last_three_bits_decimal = bin2dec(dec2bin(last_three_bits));
                
                % Subtract the last three bits from the pixel value
                temp_pixel_value_without_last_three_bits = temp_pixel_value - last_three_bits_decimal;
                
                % Add the decimal_value to the pixel value without its last three bits
                new_pixel_value = temp_pixel_value_without_last_three_bits + decimal_value;
                
                % Assign the new pixel value to smoothed_input_image(i, j)
                smoothed_input_image(i, j) = new_pixel_value;
                
              
        end

        %Do the same for it's partner.
        %disp("Partner");
       % disp("X-prime-prime for partner");
      %  disp(X_prime_prime(row+16,col));
        [row_num, col_num] = blockNumberToRowCol(X_prime_prime(row+16, col), block_size, image_size);
       
       % disp("row, col of input_image_matrix for partner");
       % disp(['Row: ', num2str(row_num)]);
       % disp(['Column: ', num2str(col_num)]);

        start_partner = 1;
        ending_partner = 3;
        for pos = 1:length(positions)
                x = positions{pos}(1);
                y = positions{pos}(2);
                i = row_num + x;
                j = col_num + y;
                % Get the value of the pixel X
                x_val = smoothed_input_image(i, j);
               
                current_bits = binary_string(start_partner:ending_partner);
                start_partner = start_partner + 3;
                ending_partner = ending_partner + 3;
              
                decimal_value = bin2dec(current_bits);
            
                % Extract the 3 least significant bits of X
                X_LSB = bitand(x_val, 7); % Equivalent to X AND (111) in binary
                
                % Calculate the difference between the watermark information and the LSBs
                v = decimal_value - X_LSB;
                
                %disp(v);
                
                % Apply the smoothing function
                if abs(v) < 5
                   smoothed_input_image(i, j) = x_val - X_LSB + v;
                elseif v >= 5
                    smoothed_input_image(i, j) = x_val - X_LSB + v - 8;
                else
                    smoothed_input_image(i, j) = x_val - X_LSB + v + 8;
                end
                
                % embed in smoothed_input_image
                 % Extract the current pixel value from smoothed_input_image
                temp_pixel_value = x_val;
                
                % Extract the last three bits of the pixel value
                last_three_bits = bitand(temp_pixel_value, 7); % Equivalent to temp_pixel_value AND 7 in binary
                
                % Convert the last three bits to decimal value
                last_three_bits_decimal = bin2dec(dec2bin(last_three_bits));
                
                % Subtract the last three bits from the pixel value
                temp_pixel_value_without_last_three_bits = temp_pixel_value - last_three_bits_decimal;
                
                % Add the decimal_value to the pixel value without its last three bits
                new_pixel_value = temp_pixel_value_without_last_three_bits + decimal_value;
                
                % Assign the new pixel value to smoothed_input_image(i, j)
                smoothed_input_image(i, j) = new_pixel_value;
                

        end
        

    end
 end
imshow(smoothed_input_image);
disp("Showing Output");
disp(smoothed_input_image);
%embedding completed
%disp(X_prime_prime);

% Given coordinates
%row = 27;
%col = 57;

% Access pixel values
%pixel_1 = input_image(row, col);
%pixel_2 = input_image(row, col + 1);
%pixel_3 = input_image(row + 1, col);
%pixel_4 = input_image(row + 1, col + 1);

% Display pixel values
%disp(['Pixel value at (27,57): ', num2str(pixel_1)]);
%disp(['Pixel value at (27,58): ', num2str(pixel_2)]);
%disp(['Pixel value at (28,57): ', num2str(pixel_3)]);
%disp(['Pixel value at (28,58): ', num2str(pixel_4)]);
% 
% row = 27;
% col = 57;
% 
% % Access pixel values
% 
% disp("After");
% pixel_1 = smoothed_input_image(row, col);
% pixel_2 = smoothed_input_image(row, col + 1);
% pixel_3 = smoothed_input_image(row + 1, col);
% pixel_4 = smoothed_input_image(row + 1, col + 1);

% % Display pixel values
% disp(['Pixel value at (27,57): ', num2str(pixel_1)]);
% disp(['Pixel value at (27,58): ', num2str(pixel_2)]);
% disp(['Pixel value at (28,57): ', num2str(pixel_3)]);
% disp(['Pixel value at (28,58): ', num2str(pixel_4)]);
% 
% 
% % Given coordinates
% row = 59;
% col = 57;
% 
% % Access pixel values
% pixel_1 = input_image(row, col);
% pixel_2 = input_image(row, col + 1);
% pixel_3 = input_image(row + 1, col);
% pixel_4 = input_image(row + 1, col + 1);
% 
% % Display pixel values
% disp(['Pixel value at (59,57): ', num2str(pixel_1)]);
% disp(['Pixel value at (59,58): ', num2str(pixel_2)]);
% disp(['Pixel value at (60,57): ', num2str(pixel_3)]);
% disp(['Pixel value at (60,58): ', num2str(pixel_4)]);
% 
% row = 59;
% col = 57;
% 
% % Access pixel values
% pixel_1 = smoothed_input_image(row, col);
% pixel_2 = smoothed_input_image(row, col + 1);
% pixel_3 = smoothed_input_image(row + 1, col);
% pixel_4 = smoothed_input_image(row + 1, col + 1);
% 
% % Display pixel values
% disp(['Pixel value at (59,57): ', num2str(pixel_1)]);
% disp(['Pixel value at (59,58): ', num2str(pixel_2)]);
% disp(['Pixel value at (60,57): ', num2str(pixel_3)]);
% disp(['Pixel value at (60,58): ', num2str(pixel_4)]);
% 
% 
% 
% %Check if the embedding algo has worked properly.
% %LET us check for a random watermark_info at watermark_matrix(10,11)
% %Corresponding X-prime-prime : (10,11), (26,11)
% disp("Watermark info 12 bits:");
% disp(watermark_matrix(10,11));
% [row_in_og_image_a,col_in_og_image_a] =blockNumberToRowCol(X_prime_prime(10,11), block_size, image_size);
% [row_in_og_image_b,col_in_og_image_b] =blockNumberToRowCol(X_prime_prime(26,11), block_size, image_size);
% 
% % Given coordinates
% row = row_in_og_image_a;
% col = col_in_og_image_a;
% 
% % Access pixel values
% pixel_1 = input_image(row, col);
% pixel_2 = input_image(row, col + 1);
% pixel_3 = input_image(row + 1, col);
% pixel_4 = input_image(row + 1, col + 1);
% 
% % Display pixel values
% disp("Before for a");
% disp(['Pixel value at top left for a: ', num2str(pixel_1)]);
% disp(['Pixel value at top right for a: ', num2str(pixel_2)]);
% disp(['Pixel value at bottom left for a: ', num2str(pixel_3)]);
% disp(['Pixel value at bottom right for a: ', num2str(pixel_4)]);
% 
% row = row_in_og_image_a;
% col = col_in_og_image_a;
% 
% % Access pixel values
% 
% disp("After for a:");
% pixel_1 = smoothed_input_image(row, col);
% pixel_2 = smoothed_input_image(row, col + 1);
% pixel_3 = smoothed_input_image(row + 1, col);
% pixel_4 = smoothed_input_image(row + 1, col + 1);
% 
% % Display pixel values
% disp(['Pixel value at top left for a): ', num2str(pixel_1)]);
% disp(['Pixel value at top right for a: ', num2str(pixel_2)]);
% disp(['Pixel value at bottom left for a: ', num2str(pixel_3)]);
% disp(['Pixel value at bottom right for a: ', num2str(pixel_4)]);
% 
% 
% % Given coordinates
% disp("For Partner: Before:");
% row = row_in_og_image_b;
% col = col_in_og_image_b;
% 
% 
% % Access pixel values
% pixel_1 = input_image(row, col);
% pixel_2 = input_image(row, col + 1);
% pixel_3 = input_image(row + 1, col);
% pixel_4 = input_image(row + 1, col + 1);
% 
% % Display pixel values
% disp(['Pixel value at top left for b: ', num2str(pixel_1)]);
% disp(['Pixel value at top right for b: ', num2str(pixel_2)]);
% disp(['Pixel value at bottom left for b: ', num2str(pixel_3)]);
% disp(['Pixel value at bottom right for b: ', num2str(pixel_4)]);
% row = row_in_og_image_b;
% col = col_in_og_image_b;
% 
% 
% % Access pixel values
% pixel_1 = smoothed_input_image(row, col);
% pixel_2 = smoothed_input_image(row, col + 1);
% pixel_3 = smoothed_input_image(row + 1, col);
% pixel_4 = smoothed_input_image(row + 1, col + 1);
% 
% % Display pixel values
% disp("After for partner");
% disp(['Pixel value at top left for b: ', num2str(pixel_1)]);
% disp(['Pixel value at top right for b: ', num2str(pixel_2)]);
% disp(['Pixel value at bottom left for b: ', num2str(pixel_3)]);
% disp(['Pixel value at bottom right for b: ', num2str(pixel_4)]);


% 2.2. Detection of tampered blocks

% Tamper detection is done in four levels. In each level we try to reduce 
% false alarms along with the localization of the tampered area. 
% We denote the image under detection as B′.

Received_image = smoothed_input_image; % Image received by the user.

valid_blocks = ones(M/2); % initialise a matrix of M/2 * M/2, Initially marked all blocks as valid.


% Level-1 detection
%For each block b, the 12-bit watermark is extracted from
%the three LSBs of each pixel. Divide the extracted watermark with the
%secret polynomial used in the embedding stage. If there is a remainder, it
%means that the block has been tampered with in the three LSBs, so mark
%b as invalid; otherwise, mark b as valid. 

%disp(watermark_matrix(1,9));
for row = 1:2:63  %1:2:31
         for col = 1:2:63  %1:2:63
             concatenated_binary_string = '';
             for pos = 1:length(positions)
                i = positions{pos}(1);
                j = positions{pos}(2);
                x = row + i;
                y = col + j;
                pixel_value = Received_image(x, y);
                pixel_binary_string = getLSBs(pixel_value);
                
                % Append the binary string to the concatenated string
                concatenated_binary_string = [concatenated_binary_string, pixel_binary_string]; %#ok<AGROW>

             end
             %Perform division with the polynomial used.
             % Perform division with the polynomial
            remainder = mod2div(concatenated_binary_string, divisor_binary);
            
            % Check if remainder is all zeros
            if all(remainder == '0')
                % Mark block as valid
                %do nothing
            else
                % Mark block as invalid
                valid_blocks((row + 1)/2, (col + 1)/2) = 0;
                disp("invalid_block detected at level 1.")
                disp(row);
                disp(col);
            end

         end
end


% Done with the Level - 1 detection.


% LEVEL 2 DETECTION.
% USING THE LOGIC OF THE RESEARCH PAPER
% For each block B which is marked valid after
% level-1 detection, check the following four triples of neighboring
% block situation: (N, NE, E), (E, SE, S), (S, SW, W), (W, NW, N). The
% 3 × 3 block-neighborhood of the central block B can be denoted by
% the compass directions:
%   NW    N     NE
%   W     B     E
%   SW    S     SE

% If all blocks in any of the four triples are marked invalid, mark block
%B invalid.

received_after_level_one = valid_blocks;

positions_top_right = {[-1, 0], [-1, 1], [0, 1]};
positions_bottom_right = { [0, 1],[1, 1], [1, 0]};
positions_bottom_left = { [1, 0],[1, -1], [0, -1]};
positions_top_left = { [0, -1],[-1, -1], [-1, 0]};






% Traverse the received_after_level_one matrix
for row = 1:M/2
    for col = 1:M/2
        if received_after_level_one(row,col) == 0
            continue;
        end
        % Check each position vector
        count_invalid = 0;
        neighbour_blocks = 0;
        for pos = 1:length(positions_top_right)
            % Get the coordinates of the neighboring block
            x = row + positions_top_right{pos}(1);
            y = col + positions_top_right{pos}(2);
            
            % Check if the coordinates are within bounds
            if x >= 1 && x <= 32 && y >= 1 && y <= 32
                neighbour_blocks = neighbour_blocks + 1;
                % Check if all neighboring blocks are invalid
                if valid_blocks(x,y) == 0
                    count_invalid = count_invalid + 1;
                   
                end  
            end
        end
        if count_invalid == neighbour_blocks && count_invalid ~= 0
            received_after_level_one(row, col) = 0;
            break;
        end
        
        count_invalid = 0;
        neighbour_blocks = 0;
        for pos = 1:length(positions_top_left)
            % Get the coordinates of the neighboring block
            x = row + positions_top_left{pos}(1);
            y = col + positions_top_left{pos}(2);
            
            % Check if the coordinates are within bounds
            if x >= 1 && x <= 32 && y >= 1 && y <= 32
                neighbour_blocks = neighbour_blocks + 1;
                % Check if all neighboring blocks are invalid
                if valid_blocks(x,y) == 0
                    count_invalid = count_invalid + 1;
                end  
            end
        end
        if count_invalid == neighbour_blocks && count_invalid ~= 0
            received_after_level_one(row, col) = 0;
            break;
        end

        count_invalid = 0;
        neighbour_blocks = 0;
        for pos = 1:length(positions_bottom_left)
            % Get the coordinates of the neighboring block
            x = row + positions_bottom_left{pos}(1);
            y = col + positions_bottom_left{pos}(2);
            
            % Check if the coordinates are within bounds
            if x >= 1 && x <= 32 && y >= 1 && y <= 32
                neighbour_blocks = neighbour_blocks + 1;
                % Check if all neighboring blocks are invalid
                if valid_blocks(x,y) == 0
                    count_invalid = count_invalid + 1;
                end  
            end
        end
        if count_invalid == neighbour_blocks && count_invalid ~= 0
            received_after_level_one(row, col) = 0;
            break;
        end

        count_invalid = 0;
        neighbour_blocks = 0;
        for pos = 1:length(positions_bottom_right)
            % Get the coordinates of the neighboring block
            x = row + positions_bottom_right{pos}(1);
            y = col + positions_bottom_right{pos}(2);
            
            % Check if the coordinates are within bounds
            if x >= 1 && x <= 32 && y >= 1 && y <= 32
                neighbour_blocks = neighbour_blocks + 1;
                % Check if all neighboring blocks are invalid
                if valid_blocks(x,y) == 0
                    count_invalid = count_invalid + 1;
                end  
            end
        end
        if count_invalid == neighbour_blocks && count_invalid ~= 0
            received_after_level_one(row, col) = 0;
            break;
        end
    end
end
% disp("Received after level one");
% 
% disp(received_after_level_one);



% Done with the Level - 2 detection.


% Level-3 detection: For each valid block after level-2 detection,
% mark this block invalid if there are five or more neighboring blocks
% in its 3 × 3 block-neighborhood that are invalid.

received_after_level_two = received_after_level_one;


for row = 1:M/2
    for col = 1:M/2
        % Check each position vector
        if received_after_level_two(row,col) == 0
            continue;
        end
        count_invalid = 0;
        neighbour_blocks = 0;
        for pos = 1:length(positions_top_right)
            % Get the coordinates of the neighboring block
            x = row + positions_top_right{pos}(1);
            y = col + positions_top_right{pos}(2);
            
            % Check if the coordinates are within bounds
            if x >= 1 && x <= 32 && y >= 1 && y <= 32
                neighbour_blocks = neighbour_blocks + 1;
                % Check if all neighboring blocks are invalid
                if received_after_level_one(x,y) == 0
                    count_invalid = count_invalid + 1;
                   
                end  
            end
        end
        if count_invalid >= 5
            received_after_level_two(row, col) = 0;
            break;
        end
        
       
        for pos = 1:length(positions_top_left)
            % Get the coordinates of the neighboring block
            x = row + positions_top_left{pos}(1);
            y = col + positions_top_left{pos}(2);
            
            % Check if the coordinates are within bounds
            if x >= 1 && x <= 32 && y >= 1 && y <= 32
                neighbour_blocks = neighbour_blocks + 1;
                % Check if all neighboring blocks are invalid
                if received_after_level_one(x,y) == 0
                    count_invalid = count_invalid + 1;
                end  
            end
        end
        if count_invalid >= 5
            received_after_level_two(row, col) = 0;
            break;
        end

        
        for pos = 1:length(positions_bottom_left)
            % Get the coordinates of the neighboring block
            x = row + positions_bottom_left{pos}(1);
            y = col + positions_bottom_left{pos}(2);
            
            % Check if the coordinates are within bounds
            if x >= 1 && x <= 32 && y >= 1 && y <= 32
                neighbour_blocks = neighbour_blocks + 1;
                % Check if all neighboring blocks are invalid
                if received_after_level_one(x,y) == 0
                    count_invalid = count_invalid + 1;
                end  
            end
        end
        if count_invalid >= 5
            received_after_level_two(row, col) = 0;
            break;
        end

       
        for pos = 1:length(positions_bottom_right)
            % Get the coordinates of the neighboring block
            x = row + positions_bottom_right{pos}(1);
            y = col + positions_bottom_right{pos}(2);
            
            % Check if the coordinates are within bounds
            if x >= 1 && x <= 32 && y >= 1 && y <= 32
                neighbour_blocks = neighbour_blocks + 1;
                % Check if all neighboring blocks are invalid
                if received_after_level_one(x,y) == 0
                    count_invalid = count_invalid + 1;
                end  
            end
        end
        if count_invalid >= 5
            received_after_level_two(row, col) = 0;
            break;
        end
    end
end


%disp(received_after_level_two);




% Done with the Level - 3 detection.
% In general, if a block has been tampered with, then there are more
% chances that the neighboring blocks have also been tampered with. Thus,
% Level 2 and Level 3 reduce the false acceptance ratio of invalid blocks and
% helps to reduce the false alarms

received_after_level_three = received_after_level_two;
% Level-4 detection: This level is to detect any changes in the five MSBs of
% each pixel of the valid blocks and to thwart VQ attacks. In this level we
% mark the block invalid only if there is a mismatch of the feature with both
% of the mapping blocks, so that false alarms can be reduced. For each valid
% block b after Level-3 detection, we set the three LSBs of each pixel in block
% b to zero. Compute the average intensity of block b, which is denoted as
% avg_b. Determine the corresponding mapping partner block x-y from the
% lookup table and perform the following task:



% If x is valid, then
% extract the 12-bit watermark from x
    % If b belongs to the {upper half of the image}, then
    % compare the five MSBs of avg b against bits 1 to 5 of the extracted
    %  watermark.
    % Else,
    %  compare the five MSBs of avg b against bits 6 to 10 of the extracted
    %  watermark.
    % If there is any mismatch in bit comparison, then mark as invalid.

    %If b is invalid, then

        %if y is valid, then
            
            %extract the 12-bit watermark from y.
            % If b belongs to the upper half of the image, then
                % compare the five MSBs of avg b _ e against bits 1 to 5 
                % of the extracted watermark.
            % Else,
                %  compare the five MSBs of avg _ b against bits 6 to 10 
                % of the extracted watermark.

        %If there is any mismatch in bit comparison, then mark as invalid.
        % Else, mark b as valid.

% Else if y is valid, then
    %extraact the 12 bit watermark from y/
    % If b belongs to the upper half, then
        %compare the fuve MSBs of avg_a against bits 1 to 5
        %of the extracted watermark.
    % Else
        %compare the five MSBs of avg_b against against bits 6 to 10 
        % of the extracted watermark
    %If there is any mismatch in bit comparison, mark b as invalid.

 
 for row = 1:32
     for col = 1:32
         if received_after_level_three(row,col) == 0
             continue;
         end
         %find the mapping blocks x, y.
           sum_pixels = 0;
         %find avg intensity, first set the LSBs zero.
          for pos = 1:length(positions)
                i = positions{pos}(1);
                j = positions{pos}(2);
                x = row + i;
                y = col + j;
                pixel_value = Received_image(x, y);
                Received_image(x,y)=setThreeLSBzero(pixel_value); % set the LSB's 0
                sum_pixels = sum_pixels + Received_image(x,y);
          end
          avg_b = sum_pixels/4;

         map_block_x = X_prime_prime(row,col);
        if row <= 16
            map_block_y = X_prime_prime(row + 16, col);
        else
            
            map_block_y = X_prime_prime(row - 16, col);
        end
        
        [start_row_x, start_col_x] = blockNumberToRowCol(map_block_x,block_size,image_size);
        [start_row_y, start_col_y] = blockNumberToRowCol(map_block_y,block_size,image_size);
        
        %map_block_x in terms of 32 x 32.
        [row_dimension_x, col_dimension_x] = blockNumberToRowColTWO(map_block_x);
         %map_block_y in terms of 32 x 32.
        [row_dimension_y, col_dimension_y] = blockNumberToRowColTWO(map_block_y);

        concatenated_binary_string_x = '';
             for pos = 1:length(positions)
                i = positions{pos}(1);
                j = positions{pos}(2);
                x = start_row_x + i;
                y = start_col_x + j;
                pixel_value = Received_image(x, y);
                pixel_binary_string = getLSBs(pixel_value);
                
                % Append the binary string to the concatenated string
                concatenated_binary_string_x = [concatenated_binary_string_x, pixel_binary_string]; %#ok<*AGROW>
             end

        concatenated_binary_string_y = '';
             for pos = 1:length(positions)
                i = positions{pos}(1);
                j = positions{pos}(2);
                x = start_row_y + i;
                y = start_col_y + j;
                pixel_value = Received_image(x, y);
                pixel_binary_string = getLSBs(pixel_value);
                
                % Append the binary string to the concatenated string
                concatenated_binary_string_y = [concatenated_binary_string_y, pixel_binary_string]; %#ok<*AGROW>
             end

        if received_after_level_three(row_dimension_x,col_dimension_x)== 1
            if row<=16
                % Extract the 5 MSBs of avg_b
                msb_avg_b = dec2bin(bitshift(avg_b, -3), 5);
                
                % Extract the first 5 bits of concatenated_binary_string_x
                bits_concatenated_binary_string_x = concatenated_binary_string_x(1:5);
                
                % Compare the MSBs of avg_b with the bits from concatenated_binary_string_x
                if ~isequal(msb_avg_b, bits_concatenated_binary_string_x)
                    received_after_level_three(row,col)=0;
                end

            else
                msb_avg_b = dec2bin(bitshift(avg_b, -3), 5);
                bits_concatenated_binary_string_x = concatenated_binary_string_x(6:10);
                 % Compare the MSBs of avg_b with the bits from concatenated_binary_string_x
                if ~isequal(msb_avg_b, bits_concatenated_binary_string_x)
                    received_after_level_three(row,col)=0;
                end

            end
        
            if received_after_level_three(row,col) == 0
                   if received_after_level_three(row_dimension_y,col_dimension_y)== 1
                        if row<=16
                            % Extract the 5 MSBs of avg_b
                            msb_avg_b = dec2bin(bitshift(avg_b, -3), 5);
                            
                            % Extract the first 5 bits of concatenated_binary_string_x
                            bits_concatenated_binary_string_y = concatenated_binary_string_y(1:5);
                            
                            % Compare the MSBs of avg_b with the bits from concatenated_binary_string_x
                            if ~isequal(msb_avg_b, bits_concatenated_binary_string_y)
                                received_after_level_three(row,col)=0;
                            else
                                received_after_level_three(row,col)=1;
                            end
            
                        else
                            msb_avg_b = dec2bin(bitshift(avg_b, -3), 5);
                            bits_concatenated_binary_string_y = concatenated_binary_string_y(6:10);
                             % Compare the MSBs of avg_b with the bits from concatenated_binary_string_x
                            if ~isequal(msb_avg_b, bits_concatenated_binary_string_y)
                                received_after_level_three(row,col)=0;
                            else
                                received_after_level_three(row,col)=1;
                            end
    
                        end
                   end
            end
        elseif received_after_level_three(row_dimension_y,col_dimension_y)== 1
            if row<=16
                % Extract the 5 MSBs of avg_b
                msb_avg_b = dec2bin(bitshift(avg_b, -3), 5);
                
                % Extract the first 5 bits of concatenated_binary_string_x
                bits_concatenated_binary_string_y = concatenated_binary_string_y(1:5);
                
                % Compare the MSBs of avg_b with the bits from concatenated_binary_string_x
                if ~isequal(msb_avg_b, bits_concatenated_binary_string_y)
                    received_after_level_three(row,col)=0;
                end

            else
                msb_avg_b = dec2bin(bitshift(avg_b, -3), 5);
                bits_concatenated_binary_string_y = concatenated_binary_string_y(6:10);
                 % Compare the MSBs of avg_b with the bits from concatenated_binary_string_x
                if ~isequal(msb_avg_b, bits_concatenated_binary_string_y)
                    received_after_level_three(row,col)=0;
               
                end

            end

        end

     end
 end
disp(received_after_level_three);







% DONE WITH ALLL FOUR DETECTIONS



% RECOVERY..

% As in Lee and Lin [20], the recovery of tampered blocks is done in two stages. Stage
% 1 is responsible for recovery from the valid mapping blocks, which is the recovery
% from the embedded watermark. In Stage 2, we recover the image by averaging the
% 3 × 3 neighborhood


for row = 1:M/2
    for col = 1:M/2
        if received_after_level_three(row,col) == 1
            continue;
        end
        x_map_block = X_prime_prime(row,col);
        if row <= 16
            map_block_y = X_prime_prime(row + 16, col);
        else
            
            map_block_y = X_prime_prime(row - 16, col);
        end
        
         candidate_block_num=0;

         [start_row_x,start_row_y]=blockNumberToRowColTWO(x_map_block);
         [start_row_x_partner,start_row_y_partner]=blockNumberToRowColTWO(map_block_y);

        if received_after_level_three(start_row_x, start_row_y) == 1 
            candidate_block_num = x_map_block;
        elseif received_after_level_three(start_row_x_partner, start_row_y_partner) == 1 
            candidate_block_num = map_block_y;
        end
        if candidate_block_num == 0
            continue;
        end
       [start_row_x,start_row_y]=blockNumberToRowCol(candidate_block_num,block_size,image_size);
       
        
         

      
     concatenated_binary_string_x = '';
     for pos = 1:length(positions)
        i = positions{pos}(1);
        j = positions{pos}(2);
        x = start_row_x + i;
        y = start_col_x + j;
        pixel_value = Received_image(x, y);
        pixel_binary_string = getLSBs(pixel_value);
        
        % Append the binary string to the concatenated string
        concatenated_binary_string_x = [concatenated_binary_string_x, pixel_binary_string]; %#ok<*AGROW>
     end
     if row <= 16
         watermark_bits = concatenated_binary_string_x(1:5);

     else
         watermark_bits = concatenated_binary_string_x(6:10);
     end
     % Pad three zeroes at the end of first_5_bits
    representative_info = [watermark_bits, '000'];
    decimal_value_to_replace = bin2dec(representative_info);
    input_image_row = (2 * row - 1);
    input_image_col = (2 * col - 1);
    for pos = 1:length(positions)
        i = positions{pos}(1);
        j = positions{pos}(2);
        x = input_image_row + i;
        y = input_image_col + j;
        smoothed_input_image(x,y) = decimal_value_to_replace;
        
     end
    received_after_level_three(row,col)=1;
        
    end
end

disp("Recovery L1 done.");
positions_two = {[-1, 0],[2,0], [0, -1], [0, 2], [-1, 1],[1,2],[2,-1],[2,1],[2,2],[-1,-1],[-1,1],[-1,2]};
for row = 1:M/2
    for col = 1:M/2
        if received_after_level_three(row,col) == 1
            continue;
        end

        smoothed_image_row = (2 * row) - 1;
        smoothed_image_col = (2 * col) - 1;
        sum_neighbours= 0;
        neighbour_count = 0;
        for pos = 1:length(positions_two)
            i = positions{pos_two}(1);
            j = positions{pos_two}(2);
            x = smoothed_image_row + i;
            y = smoothed_image_row + j;
            if x>=0 && x<=64 && y>=0 && y<=64
                sum_neighbours = smoothed_input_image(x,y);
                neighbour_count = neighbour_count + 1;
            end
        end
        avg_neighbour = sum_neighbours/neighbour_count;
        for pos = 1:length(positions_two)
            i = positions{pos_two}(1);
            j = positions{pos_two}(2);
            x = smoothed_image_row + i;
            y = smoothed_image_row + j;
            smoothed_input_image(x,y) = avg_neighbour;
                
            
        end
        received_after_level_three = 1;
    end

 end

 disp("Recovery after L2 done");


function result = xor1(a, b)
    % Initialize result
    result = '';
    
    % Traverse all bits, if bits are same, then XOR is 0, else 1
    for i = 1:length(b)
        if a(i) == b(i)
            result = strcat(result, '0');
        else
            result = strcat(result, '1');
        end
    end
end

function remainder = mod2div(dividend, divisor)
    % Number of bits to be XORed at a time.
    dividend = dividend(find(dividend == '1', 1):end);
    if length(dividend) < length(divisor)
        remainder = '00'; % Return '00' as the remainder
        return;
    end
    pick = length(divisor);
    % Slicing the dividend to the appropriate length for a particular step
    tmp = dividend(1:pick);
    n = length(dividend);
    while pick < n   
            s1 = xor1(divisor, tmp);
            tmp = s1;
            tmp = tmp(find(tmp == '1', 1):end);  % Remove leading zeros
            % Increment pick to move further
            while length(tmp) < length(divisor)
                % Concatenate dividend bits until length matches divisor's length
                pick = pick + 1;   
                if pick <= length(dividend)
                    % Concatenate dividend bits until length matches divisor's length
                    tmp = strcat(tmp, dividend(pick));
                else
                    break; % Exit the loop if pick exceeds the length of the dividend
                end
            end     
    end
    if(pick<=length(dividend))
        tmp = tmp(find(tmp == '1', 1):end);
                s1 = xor1(divisor, tmp);
                tmp = s1;
    
    end
    if length(tmp) < 5
        tmp = [repmat('0', 1, 5 - length(tmp)) tmp];
    end
      remainder = tmp;
end


function [row, col] = blockNumberToRowCol(block_number, block_size, image_size)
    % Calculate the number of blocks in each dimension
    num_blocks_row = image_size(1) / block_size(1);
    num_blocks_col = image_size(2) / block_size(2);

    % Calculate the row and column indices of the block
    col = mod(block_number - 1, num_blocks_col) * block_size(2) + 1;
    row = ceil(block_number / num_blocks_row);
    row = 2*row - 1;
end

function lsb_chars = getLSBs(pixel_value)
    % Extract the 3 LSBs of the pixel value
    lsb_bits = bitand(pixel_value, 7);  % Equivalent to pixel_value AND 111 in binary
    
    % Convert the LSBs to binary string format
    lsb_chars = dec2bin(lsb_bits, 3);
end


function modified = setThreeLSBzero(pixel_value)
     last_three_bits = bitand(pixel_value, 7); % Equivalent to temp_pixel_value AND 7 in binary
                
    % Convert the last three bits to decimal value
    last_three_bits_decimal = bin2dec(dec2bin(last_three_bits));
    
    % Subtract the last three bits from the pixel value
    modified = pixel_value - last_three_bits_decimal;
    
end


function [row, col] = blockNumberToRowColTWO(block_number)
    % Calculate the row and column indices of the block
    col = mod(block_number, 32);
    if col == 0
        col = 32;
    end
    row = ceil(block_number /32);
    
end

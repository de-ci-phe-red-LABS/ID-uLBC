int: RD;

% 32 nibble permutation for uLBC
array[0..31] of int: post_perm = array1d(0..31, [0, 25, 18, 15, 4, 1, 30, 19, 8, 29, 22, 7, 12, 9, 26, 23, 20, 5, 10, 31, 16, 13, 2, 27, 28, 17, 6, 3, 24, 21, 14, 11]);

% #############################################################################################################################################
% Constraints for EU (Upper Trail)
array[0..RD, 0..31] of var 0..3: AXU;
array[0..RD, 0..31] of var -2..15: DXU;
constraint forall(i in 0..RD, j in 0..31) (
    link_difference(AXU[i,j], DXU[i,j])
);

array[0..(RD - 1), 0..31] of var 0..3: AYU;
array[0..(RD - 1), 0..31] of var -2..15: DYU;
constraint forall(i in 0..(RD - 1), j in 0..31) (
    link_difference(AYU[i,j], DYU[i,j])
);

array[0..(RD - 1), 0..31] of var 0..3: AZU;
array[0..(RD - 1), 0..31] of var -2..15: DZU;
constraint forall(i in 0..(RD - 1), j in 0..31) (
    link_difference(AZU[i,j], DZU[i,j])
);

% exclude all-zero input difference
constraint sum(i in 0..31)(AXU[0, i]) != 0;

% S-box
constraint forall(r in 0..(RD - 1), i in 0..31)
(
    sbox(AXU[r, i], AYU[r, i])
);

% PostPerm (Linear Layer - Permutation)
constraint forall(r in 0..(RD - 1), i in 0..31)
(
    AYU[r, post_perm[i]] = AZU[r, i] /\
    DYU[r, post_perm[i]] = DZU[r, i]
);

% MixColumn (Linear Layer - Matrix)
constraint forall(r in 0..(RD - 1), i in 0..7)
(
    mix_column(AZU[r, 4*i + 0], DZU[r, 4*i + 0],
               AZU[r, 4*i + 1], DZU[r, 4*i + 1],
               AZU[r, 4*i + 2], DZU[r, 4*i + 2],
               AZU[r, 4*i + 3], DZU[r, 4*i + 3],
               AXU[r + 1, 4*i + 0], DXU[r + 1, 4*i + 0],
               AXU[r + 1, 4*i + 1], DXU[r + 1, 4*i + 1],
               AXU[r + 1, 4*i + 2], DXU[r + 1, 4*i + 2],
               AXU[r + 1, 4*i + 3], DXU[r + 1, 4*i + 3])
);

% #############################################################################################################################################
% Constraints for EL (Lower Trail)

array[0..RD, 0..31] of var 0..3: AXL;
array[0..RD, 0..31] of var -2..15: DXL;
constraint forall(i in 0..RD, j in 0..31) (
    link_difference(AXL[i,j], DXL[i,j])
);

array[0..(RD - 1), 0..31] of var 0..3: AYL;
array[0..(RD - 1), 0..31] of var -2..15: DYL;
constraint forall(i in 0..(RD - 1), j in 0..31) (
    link_difference(AYL[i,j], DYL[i,j])
);

array[0..(RD - 1), 0..31] of var 0..3: AZL;
array[0..(RD - 1), 0..31] of var -2..15: DZL;
constraint forall(i in 0..(RD - 1), j in 0..31) (
    link_difference(AZL[i, j], DZL[i, j])
);

% exclude all-zero output difference
constraint sum(i in 0..31)(AXL[RD, i]) != 0;

% MixColumn (Backward)
constraint forall(r in 0..(RD - 1), i in 0..7)
(
    mix_column(AXL[r + 1, 4*i + 0], DXL[r + 1, 4*i + 0],
               AXL[r + 1, 4*i + 1], DXL[r + 1, 4*i + 1],
               AXL[r + 1, 4*i + 2], DXL[r + 1, 4*i + 2],
               AXL[r + 1, 4*i + 3], DXL[r + 1, 4*i + 3],
               AZL[r, 4*i + 0], DZL[r, 4*i + 0],
               AZL[r, 4*i + 1], DZL[r, 4*i + 1],
               AZL[r, 4*i + 2], DZL[r, 4*i + 2],
               AZL[r, 4*i + 3], DZL[r, 4*i + 3])
);

% PostPerm (Backward)
constraint forall(r in 0..(RD - 1), i in 0..31)
(
    AZL[r, i] = AYL[r, post_perm[i]] /\
    DZL[r, i] = DYL[r, post_perm[i]]
);

% S-box (backward)
constraint forall(r in 0..(RD - 1), i in 0..31)
(
    sbox(AYL[r, i], AXL[r, i])
);

% #############################################################################################################################################
% Contradictions

array[0..RD] of var 0..32: contradict1;
array[0..RD] of var 0..32: contradict2;
constraint forall(r in 0..RD)
(
    contradict1[r] = sum(i in 0..31)(bool2int((AXU[r, i] + AXL[r, i] > 0) /\ (AXU[r, i] + AXL[r, i] < 3) /\ (AXU[r, i] != AXL[r, i]))) /\
    contradict2[r] = sum(i in 0..31)(bool2int(AXU[r, i] == 1 /\ AXL[r, i] == 1 /\ (DXU[r, i] != DXL[r, i])))
);

array[0..(RD - 1)] of var 0..32: contradict3;
array[0..(RD - 1)] of var 0..32: contradict4;
constraint forall(r in 0..(RD - 1))
(
    contradict3[r] = sum(i in 0..31)(bool2int((AYU[r, i] + AYL[r, i] > 0) /\ (AYU[r, i] + AYL[r, i] < 3) /\ (AYU[r, i] != AYL[r, i]))) /\
    contradict4[r] = sum(i in 0..31)(bool2int(AYU[r, i] == 1 /\ AYL[r, i] == 1 /\ (DYU[r, i] != DYL[r, i])))
);
var 0..1000: contradict;
constraint contradict = sum(r in 0..(RD - 1))(contradict1[r] + contradict2[r] + contradict3[r] + contradict4[r]) + contradict1[RD] + contradict2[RD];
constraint contradict > 0;

var int: active_bytes = sum(i in 0..31)(bool2int(AXU[0, i] > 0) + bool2int(AXL[RD, i] > 0));
solve minimize active_bytes;

% #############################################################################################################################################
% Auxiliary Functions

predicate link_difference(var 0..3: diff_pattern, var -2..15: diff_value) =  
    if (diff_pattern == 0) then diff_value == 0
    elseif (diff_pattern == 1) then diff_value > 0
    elseif (diff_pattern == 2) then diff_value == -1
    else diff_value = -2
    endif
;

predicate sbox(var 0..3: diff_in, var 0..3: diff_out) = 
    diff_out != 1 /\
    (diff_in + diff_out) in {0, 3, 4, 6} /\
    diff_out >= diff_in /\
    (diff_out - diff_in) <= 1
;

predicate xor_operation(var 0..3: diff_a, var -2..15: dvalue_a, 
                        var 0..3: diff_b, var -2..15: dvalue_b, 
                        var 0..3: diff_c, var -2..15: dvalue_c) = 
    if (diff_a + diff_b > 2) then 
        (diff_c = 3) /\ (dvalue_c = -2)
    elseif (diff_a + diff_b = 1) then
        (diff_c = 1) /\ (dvalue_c = dvalue_a + dvalue_b)
    elseif ((diff_a == 0) /\ (diff_b = 0)) then
        (diff_c = 0) /\ (dvalue_c = 0)
    elseif (dvalue_a + dvalue_b < 0) then
        (diff_c = 2) /\ (dvalue_c = -1)
    elseif (dvalue_a == dvalue_b) then
        (diff_c = 0) /\ (dvalue_c = 0)
    else
        (diff_c = 1) /\ 
        dvalue_c = (bool2int((dvalue_a mod 2) != (dvalue_b mod 2)) + 
        bool2int(((dvalue_a div 2) mod 2) != ((dvalue_b div 2) mod 2)) * 2 + 
        bool2int(((dvalue_a div 4) mod 2) != ((dvalue_b div 4) mod 2)) * 4 +
        bool2int(((dvalue_a div 8) mod 2) != ((dvalue_b div 8) mod 2)) * 8)
    endif
;

predicate mix_column(var 0..3: diff_in0, var -2..15: dvalue_in0,
                     var 0..3: diff_in1, var -2..15: dvalue_in1,
                     var 0..3: diff_in2, var -2..15: dvalue_in2,
                     var 0..3: diff_in3, var -2..15: dvalue_in3,
                     var 0..3: diff_out0, var -2..15: dvalue_out0,
                     var 0..3: diff_out1, var -2..15: dvalue_out1,
                     var 0..3: diff_out2, var -2..15: dvalue_out2,
                     var 0..3: diff_out3, var -2..15: dvalue_out3) =
    let {
        var 0..3: tmp_12; var -2..15: dtmp_12;
        var 0..3: tmp_02; var -2..15: dtmp_02;
        var 0..3: tmp_01; var -2..15: dtmp_01;
    } in (
        % y0 = x1 ^ x2 ^ x3
        xor_operation(diff_in1, dvalue_in1, diff_in2, dvalue_in2, tmp_12, dtmp_12) /\
        xor_operation(tmp_12, dtmp_12, diff_in3, dvalue_in3, diff_out0, dvalue_out0) /\
        
        % y1 = x0 ^ x2 ^ x3
        xor_operation(diff_in0, dvalue_in0, diff_in2, dvalue_in2, tmp_02, dtmp_02) /\
        xor_operation(tmp_02, dtmp_02, diff_in3, dvalue_in3, diff_out1, dvalue_out1) /\
        
        % y2 = x0 ^ x1 ^ x3
        xor_operation(diff_in0, dvalue_in0, diff_in1, dvalue_in1, tmp_01, dtmp_01) /\
        xor_operation(tmp_01, dtmp_01, diff_in3, dvalue_in3, diff_out2, dvalue_out2) /\
        
        % y3 = x0 ^ x1 ^ x2
        xor_operation(tmp_01, dtmp_01, diff_in2, dvalue_in2, diff_out3, dvalue_out3)
    )
;

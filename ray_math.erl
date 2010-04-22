-module(ray_math).
-compile(export_all).

%%%% Vector3 is defined as a 3-element tuple: {X, Y, Z}.
%%%% Color4 is defined as a 4-element tuple: {R, G, B, A}.
%%%% Matrix44 is defined as a 16-element tuple, as follows:
%%%% {	A00, A01, A02, A03,
%%%%	A10, A11, A12, A13,
%%%%	A20, A21, A22, A23,
%%%%	A30, A31, A32, A33	}

%%% General functions

%% large number - what we'll consider as infinity
large_number() -> 1.0e10.

%% signum - 1 if the number is positive, -1 if negative, 0 if equals 0.
sign(Number) when Number > 0 -> 1;
sign(Number) when Number < 0 -> -1;
sign(Number) when Number == 0 -> 0.

%% floating-point division with correct-signed 'infinities' instead of division-by-zero error
%% TODO compare with the original source of this idea to see if there is anyting missing 
fdiv(Numerator, Denominator) when Denominator == 0 -> sign(Numerator) * large_number();
fdiv(Numerator, Denominator) -> Numerator / Denominator.

%% Linear interpolation for real numbers. Scale == 0 means A, Scale == 1 means B.
%% Values out of 0..1 range lead to linear extrapolation ;).
lerp(A, B, Scale) -> (1-Scale)*A + Scale*B.

%%% Polynomial equation solvers
%% 1'st degree: ax + b = 0
%% returns the root or an atom on error
linear_equation_root(A, _) when A == 0 -> 'This is not a linear equation.';
linear_equation_root(A, B) -> fdiv(-B, A).

%% 2'nd degree: ax^2 + bx + c = 0
%% return list of real roots (may be empty!) or an atom on error
quadratic_equation_root(A, B, _) when (A == 0) and (B == 0) -> 'This is not a quadratic equation';
quadratic_equation_root(A, B, C) when A == 0 -> [linear_equation_root(B, C)];	%FIXME maybe an error should be signalized here? 
quadratic_equation_root(A, B, C) ->
	case (B*B - 4*A*C) of
		Delta when Delta < 0 -> [];
		Delta when Delta == 0 -> [fdiv(-B, 2*A)];
		Delta when Delta > 0 -> [fdiv(-B + math:sqrt(Delta), 2*A), fdiv(-B - math:sqrt(Delta), 2*A)]
	end.

%%% Vector3 functions
%% A + B
vector3_add( {Ax, Ay, Az}, {Bx, By, Bz} ) -> {Ax + Bx, Ay + By, Az + Bz}.

%% A - B.
vector3_diff( {Ax, Ay, Az}, {Bx, By, Bz} ) -> {Ax - Bx, Ay - By, Az - Bz}.

%% k*A
vector3_scale( K, {Ax, Ay, Az} ) -> { K*Ax, K*Ay, K*Az}.

%% (1/k)*A
vector3_scale_inv( K, {Ax, Ay, Az} ) -> vector3_scale(fdiv(1, K), {Ax, Ay, Az}).

%% square of vector's length
vector3_length_sqr( {Ax, Ay, Az} ) -> Ax*Ax + Ay*Ay + Az*Az.

%% vector length
vector3_length( {Ax, Ay, Az} ) -> math:sqrt(vector3_length_sqr( {Ax, Ay, Az} )).

%% normalize - reduce to unit length
vector3_normalize( {Ax, Ay, Az} ) ->
	vector3_scale( fdiv(1, vector3_length( {Ax, Ay, Az} )), {Ax, Ay, Az} ).
	
%% dot product of two vectors
vector3_dot_product( {Ax, Ay, Az}, {Bx, By, Bz} ) -> Ax*Bx + Ay*By + Az*Bz.

%% cross product of two vectors
vector3_cross_product( {Ax, Ay, Az}, {Bx, By, Bz} ) ->
	{	(Ay * Bz) - (Az * By),
		(Az * Bx) - (Ax * Bz),
		(Ax * By) - (Ay * Bx)	}.

%% reflection of vector A from surface with normal vector
%% pointing in the same direction as vector B
vector3_reflection(A, B) ->
	vector3_reflection_normalized(A, vector3_normalize(B)).

%% reflection of vector A from surface with normal vector N
vector3_reflection_normalized(A, N) -> 
	vector3_scale(vector3_length(A), vector3_diff(A, vector3_scale(2.0 * vector3_dot_product(A, N), N))).
	
%% Linear interpolation for 3-D vectors. Scale == 0 means A, Scale == 1 means B.
%% Values out of 0..1 range lead to linear extrapolation ;).
vector3_lerp(A, B, Scale) -> vector3_add(vector3_scale(1 - Scale, A), vector3_scale(Scale, B)).
	
%%% Color4 functions
%% A + B
color4_add( {Ar, Ag, Ab, Aa}, {Br, Bg, Bb, Ba} ) -> {Ar + Br, Ag + Bg, Ab + Bb, Aa + Ba}.

%% A - B.
color4_diff( {Ar, Ag, Ab, Aa}, {Br, Bg, Bb, Ba} ) -> {Ar - Br, Ag - Bg, Ab - Bb, Aa - Ba}.

%% k*A
color4_scale( K, {Ar, Ag, Ab, Aa} ) -> { K*Ar, K*Ag, K*Ab, K*Aa}.

%% (1/k)*A
color4_scale_inv( K, Color ) -> color4_scale(fdiv(1, K), Color).

%% Linear interpolation for RGBA colors. Scale == 0 means A, Scale == 1 means B.
%% Values out of 0..1 range lead to linear extrapolation ;).
color4_lerp(A, B, Scale) -> color4_add(color4_scale(1 - Scale, A), color4_scale(Scale, B)).

%% Convert (i, j) matrix indices to the number of tuple element.
matrix4_element(I, J) -> I*4+J.

%% Get identity matrix.
matrix4_identity() ->
		{	1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1	}.

%% A + B
matrix4_add( {	A00, A01, A02, A03,
				A10, A11, A12, A13,
				A20, A21, A22, A23,
				A30, A31, A32, A33	},
			{	B00, B01, B02, B03,
				B10, B11, B12, B13,
				B20, B21, B22, B23,
				B30, B31, B32, B33	} ) ->
	{	A00 + B00, A01 + B01, A02 + B02, A03 + B03,
		A10 + B10, A11 + B11, A12 + B12, A13 + B13,
		A20 + B20, A21 + B21, A22 + B22, A23 + B23,
		A30 + B30, A31 + B31, A32 + B32, A33 + B33	}.
		
%% A - B. 
matrix4_diff( {	A00, A01, A02, A03,
				A10, A11, A12, A13,
				A20, A21, A22, A23,
				A30, A31, A32, A33	},
			{	B00, B01, B02, B03,
				B10, B11, B12, B13,
				B20, B21, B22, B23,
				B30, B31, B32, B33	} ) ->
	{	A00 - B00, A01 - B01, A02 - B02, A03 - B03,
		A10 - B10, A11 - B11, A12 - B12, A13 - B13,
		A20 - B20, A21 - B21, A22 - B22, A23 - B23,
		A30 - B30, A31 - B31, A32 - B32, A33 - B33	}.

%% A*B
%% This function is rewritten from C++ code of Regedit's CommonLib - http://regedit.gamedev.pl		
matrix4_multiply( {	A00, A01, A02, A03,
				A10, A11, A12, A13,
				A20, A21, A22, A23,
				A30, A31, A32, A33	},
			{	B00, B01, B02, B03,
				B10, B11, B12, B13,
				B20, B21, B22, B23,
				B30, B31, B32, B33	} ) ->


	{	A00 * B00 + A01 * B10 + A02 * B20 + A03 * B30,
		A00 * B01 + A01 * B11 + A02 * B21 + A03 * B31,
		A00 * B02 + A01 * B12 + A02 * B22 + A03 * B32,
		A00 * B03 + A01 * B13 + A02 * B23 + A03 * B33,
		
		A10 * B00 + A11 * B10 + A12 * B20 + A13 * B30,
		A10 * B01 + A11 * B11 + A12 * B21 + A13 * B31,
		A10 * B02 + A11 * B12 + A12 * B22 + A13 * B32,
		A10 * B03 + A11 * B13 + A12 * B23 + A13 * B33,
		
		A20 * B00 + A21 * B10 + A22 * B20 + A23 * B30,
		A20 * B01 + A21 * B11 + A22 * B21 + A23 * B31,
		A20 * B02 + A21 * B12 + A22 * B22 + A23 * B32,
		A20 * B03 + A21 * B13 + A22 * B23 + A23 * B33,
		
		A30 * B00 + A31 * B10 + A32 * B20 + A33 * B30,
		A30 * B01 + A31 * B11 + A32 * B21 + A33 * B31,
		A30 * B02 + A31 * B12 + A32 * B22 + A33 * B32,
		A30 * B03 + A31 * B13 + A32 * B23 + A33 * B33	}.
	
%% det(A), aka. |A|.
%% This function is rewritten from C++ code of Regedit's CommonLib - http://regedit.gamedev.pl
matrix4_det( {	A00, A01, A02, A03,
				A10, A11, A12, A13,
				A20, A21, A22, A23,
				A30, A31, A32, A33	} ) ->

	(A00 * A11 - A10 * A01) * (A22 * A33 - A32 * A23) - (A00 * A21 - A20 * A01) * (A12 * A33 - A32 * A13) +
	(A00 * A31 - A30 * A01) * (A12 * A23 - A22 * A13) + (A10 * A21 - A20 * A11) * (A02 * A33 - A32 * A03) -
	(A10 * A31 - A30 * A11) * (A02 * A23 - A22 * A03) + (A20 * A31 - A30 * A21) * (A02 * A13 - A12 * A03).

%% Invert matrix.
%% This function is rewritten from C++ code of Regedit's CommonLib - http://regedit.gamedev.pl		
matrix4_inv( {	A00, A01, A02, A03,
				A10, A11, A12, A13,
				A20, A21, A22, A23,
				A30, A31, A32, A33	} ) ->
	D = fdiv(1.0, matrix4_det( {	A00, A01, A02, A03,
				A10, A11, A12, A13,
				A20, A21, A22, A23,
				A30, A31, A32, A33	} )), 
	{	D * (A11 * (A22 * A33 - A32 * A23) + A21 * (A32 * A13 - A12 * A33) + A31 * (A12 * A23 - A22 * A13)),
		D * (A21 * (A02 * A33 - A32 * A03) + A31 * (A22 * A03 - A02 * A23) + A01 * (A32 * A23 - A22 * A33)),
		D * (A31 * (A02 * A13 - A12 * A03) + A01 * (A12 * A33 - A32 * A13) + A11 * (A32 * A03 - A02 * A33)),
		D * (A01 * (A22 * A13 - A12 * A23) + A11 * (A02 * A23 - A22 * A03) + A21 * (A12 * A03 - A02 * A13)),
		D * (A12 * (A20 * A33 - A30 * A23) + A22 * (A30 * A13 - A10 * A33) + A32 * (A10 * A23 - A20 * A13)),
		D * (A22 * (A00 * A33 - A30 * A03) + A32 * (A20 * A03 - A00 * A23) + A02 * (A30 * A23 - A20 * A33)),
		D * (A32 * (A00 * A13 - A10 * A03) + A02 * (A10 * A33 - A30 * A13) + A12 * (A30 * A03 - A00 * A33)),
		D * (A02 * (A20 * A13 - A10 * A23) + A12 * (A00 * A23 - A20 * A03) + A22 * (A10 * A03 - A00 * A13)),
		D * (A13 * (A20 * A31 - A30 * A21) + A23 * (A30 * A11 - A10 * A31) + A33 * (A10 * A21 - A20 * A11)),
		D * (A23 * (A00 * A31 - A30 * A01) + A33 * (A20 * A01 - A00 * A21) + A03 * (A30 * A21 - A20 * A31)),
		D * (A33 * (A00 * A11 - A10 * A01) + A03 * (A10 * A31 - A30 * A11) + A13 * (A30 * A01 - A00 * A31)),
		D * (A03 * (A20 * A11 - A10 * A21) + A13 * (A00 * A21 - A20 * A01) + A23 * (A10 * A01 - A00 * A11)),
		D * (A10 * (A31 * A22 - A21 * A32) + A20 * (A11 * A32 - A31 * A12) + A30 * (A21 * A12 - A11 * A22)),
		D * (A20 * (A31 * A02 - A01 * A32) + A30 * (A01 * A22 - A21 * A02) + A00 * (A21 * A32 - A31 * A22)),
		D * (A30 * (A11 * A02 - A01 * A12) + A00 * (A31 * A12 - A11 * A32) + A10 * (A01 * A32 - A31 * A02)),
		D * (A00 * (A11 * A22 - A21 * A12) + A10 * (A21 * A02 - A01 * A22) + A20 * (A01 * A12 - A11 * A02))	}.

%% Transform a vector by a matrix, with translation
matrix4_transform_vector3_point( {	A00, A01, A02, _,
									A10, A11, A12, _,
									A20, A21, A22, _,
									A30, A31, A32, _	},
								{ Vx, Vy, Vz } ) ->
	{	Vx * A00 + Vy * A10 + Vz * A20 + A30,
		Vx * A01 + Vy * A11 + Vz * A21 + A31,
		Vx * A02 + Vy * A12 + Vz * A22 + A32	}.

%% Transform a vector by a matrix, without translation		
matrix4_transform_vector3_direction ( {	A00, A01, A02, _,
										A10, A11, A12, _,
										A20, A21, A22, _,
										_, _, _, _	},
									{ Vx, Vy, Vz } ) ->
	{	Vx * A00 + Vy * A10 + Vz * A20,
		Vx * A01 + Vy * A11 + Vz * A21,
		Vx * A02 + Vy * A12 + Vz * A22	}.

%% Generate a translation matrix from a vector.		
matrix4_translation( {X, Y, Z} ) ->
	{	1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		X, Y, Z, 1	}.

%% Generate a scaling matrix from a vector.		
matrix4_scale( { X, Y, Z} ) ->
	{	X, 0, 0, 0,
		0, Y, 0, 0,
		0, 0, Z, 0,
		0, 0, 0, 1	}.

%% Generate a rotation matrix for rotation around X axis with given angle (in radians).		
matrix4_rotationX(Angle) ->
	Sin = math:sin(Angle),
	Cos = math:cos(Angle),
	
	{	1, 0, 0, 0,
		0, Cos, Sin, 0,
		0, -Sin, Cos, 0,
		0, 0, 0, 1	}.
		
%% Generate a rotation matrix for rotation around Y axis with given angle (in radians).		
matrix4_rotationY(Angle) ->
	Sin = math:sin(Angle),
	Cos = math:cos(Angle),
	
	{	Cos, 0, -Sin, 0,
		0, 1, 0, 0,
		Sin, 0, Cos, 0,
		0, 0, 0, 1	}.
		
%% Generate a rotation matrix for rotation around Z axis with given angle (in radians).		
matrix4_rotationZ(Angle) ->
	Sin = math:sin(Angle),
	Cos = math:cos(Angle),
	
	{	Cos, Sin, 0, 0,
		-Sin, Cos, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1	}.
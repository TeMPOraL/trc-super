-module(raytracer).

%% real (public) exports:
-export([create_viewport/1, save_rendering_to_24bpp_bitmap/3]).

%% temporary (debug) exports:
-export([get_kwlist_param/2, bitmap_get_bmp_header/2,
			bitmap_encode_rgb/1, bitmap_write_file/3, bitmap_encode_rgb_list/1,
			bitmap_encode_rgb_list_tail/2, do_raytracing_tail/3, do_raytracing/1]).


%% Get a value from key in keyword list
get_kwlist_param([{Name, Value} | _], Name) -> Value;
get_kwlist_param([_ | Rest], Name) -> get_kwlist_param(Rest, Name);
get_kwlist_param(_,_) -> 'Error - parameter not found'.

%% rendering outputs
%% render_to_bitmap(OutW, OutH, OutData) -> ...

%% Following code implements saving to a 24bit color bitmap file.

% Generate binary BMP file header, given image width and height.
bitmap_get_bmp_header(Width, Height) ->
					<<	"BM",		% magic
						0:32/little,		% size (dword)
						0:16/little,		% reserved 1 (word)
						0:16/little,		% reserved 2 (word)
						16#36:32/little,	% offset to bitmap data (dword)
						16#28:32/little,	% size of something (dword)
						Width:32/little,	% width (in pixels) (dword)
						Height:32/little,	% height (in pixels) (dword)
						1:16/little,		% planes in target device (word)
						16#18:16/little,	% bits per pixel (word)
						0:(6*32)/little		% some other stuff (6 dwords)
					>>.

% Encode RGB color in binary BMP format.					
bitmap_encode_rgb( {Red,Green,Blue} ) ->
					<< Blue:8, Green:8, Red:8 >>.
					
% Helper function for bitmap_encode_rgb_list - does tail recursion					
bitmap_encode_rgb_list_tail(ListElement, [Head | Rest]) ->
	[bitmap_encode_rgb_list_tail(Head, Rest) | bitmap_encode_rgb(ListElement)];	% reverse order (BMP)
bitmap_encode_rgb_list_tail(ListElement, []) -> [bitmap_encode_rgb(ListElement)].
										
% Use this function to encode list of pixels into BMP binary format.
bitmap_encode_rgb_list([Head | Rest]) -> list_to_binary(bitmap_encode_rgb_list_tail(Head, Rest)).
					
bitmap_write_file(FileName, Header, Content) ->
	case file:open(FileName, [write, binary, raw]) of
		{ok, Handle} ->
			file:write(Handle, Header),
			file:write(Handle, Content),
			file:close(Handle),
			{ok};
		{error, Reason} ->
			{error, 'File could not be opened', Reason}	
	end.

%% scene object functions

%% ray_intersect(Ray, Object) -> {Object, Coord}...
%% ray_intersection_material_params(Ray, Object) -> Color, recursive

%% cast_ray(Direction) -> Color, recursive

do_raytracing_tail(VWidth, VHeight, 1) -> [{1, 2, 3}];
do_raytracing_tail(VWidth, VHeight, N) when (N > 0) -> [{N, 2*N, 5*N} | do_raytracing_tail(VWidth, VHeight, N-1)];
do_raytracing_tail(VWidth, VHeight, _) -> 'Error - it seems that N has gone negative.'.


do_raytracing([VWidth, VHeight | Rest]) -> do_raytracing_tail(VWidth, VHeight, VWidth*VHeight).


%%% Public functions
%% create_viewport/1 - create a rendering viewport. Input argument should be a
%% list of pairs: {parameter name, parameter value}. Accepted parameters are:
%% width - Viewport width in pixels.
%% height - Viewport height in pixels.
%% ? color - "yes" or "no", use color in rendering.
%% fov - field of view in angles. WARNING, this parameter will overwrite distance.
%% distance - rendering distance in units. WARNING, this parameter will overwrite
%%			  the FOV
create_viewport([]) -> 'Error - no viewport parameters!';
create_viewport(LIST) ->
	[ get_kwlist_param(LIST, width), get_kwlist_param(LIST, height), get_kwlist_param(LIST, mode) ].

%% save_rendering_to_24bpp_bitmap/3 - save rendering data to a 24-bit color bitmap file
%% Name - file name
%% VWidth, VHeight - viewport width/height, used to render Data
%% Data - list of colors - rendering result
%%
%% FIXME - this function seems to write pixels in the reverse order, somehow.
%% NOTE - this function accepts only VWidth / VHeight that are dividable by 4. It'll probably be fixed one day. 
save_rendering_to_24bpp_bitmap(Name, [VWidth, VHeight | _ ], Data) when ((VWidth rem 4) == 0) and ((VHeight rem 4) == 0) ->
	bitmap_write_file(Name, bitmap_get_bmp_header(VWidth, VHeight), bitmap_encode_rgb_list(Data)).


%% scene management
%% sphere
%% cube
%% ect
%% add_object(Scene, Object) -> ...

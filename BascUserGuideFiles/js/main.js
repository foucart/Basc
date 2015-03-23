
var initialized = false;

$(document).ready(function() {
    
    if(!initialized){
        
        //hide the interactive guide "hide example" button
        $("#hide_button_div").slideUp(0);

        //prepare interactive guide
        $("#ig_polyn_polynom_notcont_without").slideUp();
        $("#ig_polyn_polynom_notcont_with").slideUp();
        $("#ig_polyn_spline_notcont_without").slideUp();
        $("#ig_polyn_spline_cont_without").slideUp();
        $("#ig_polyn_spline_diff_without").slideUp();
        $("#ig_polyn_spline_cont_with").slideUp();
        $("#ig_polyn_spline_diff_with").slideUp();
        $("#ig_spline_polynom_notcont_without").slideUp();
        $("#ig_spline_polynom_notcont_with").slideUp();
        $("#ig_spline_spline_notcont_without").slideUp();
        $("#ig_spline_spline_notcont_with").slideUp();
        $("#ig_spline_spline_cont_without").slideUp();
        $("#ig_spline_spline_diff_without").slideUp();
        $("#ig_spline_spline_cont_with").slideUp();
        $("#ig_spline_spline_diff_with").slideUp();
        $("#ig_polyn_spline_notcont_with").slideUp();

        //get current size of elements


        //bind an event handler for all the visibility toggling buttons
    //    var body_content_container_array = $(document).find('div.body_content_container');
    //    for (i = 0; i < body_content_container_array.length; i++) {
    //        var vis_toggle_elements = body_content_container_array[i].find('div.vis_toggle');
    //        for (j = 0; j < vis_toggle_elements.length; j++) {
    //            vis_toggle_elements[j].click = function(){
    //                //vis_toggle_elements.parentElement
    //                vis_toggle_elements[j].toggleClass('rot180');
    //            };
    //        }
    //    }
        $(document).find('div.body_content_container_container').find('div.content_content_box').slideUp();

        $('#examples_norm').slideDown();
        $('#examples_weight').slideDown();
        $('#examples_breakpoints').slideDown();
        $('#examples_smoothness').slideDown();
        $('#examples_parity').slideDown();
        $('#examples_interpolation0').slideDown();
        $('#examples_interpolation1').slideDown();
        $('#examples_upper_range').slideDown();
        $('#examples_lower_range').slideDown();
        $('#examples_shape').slideDown();
        $('#examples_obj_fun').slideDown();
        $('#examples_side').slideDown();
        $('#examples_sign').slideDown();
        $('#examples_monotonicity').slideDown();
        $('#examples_concavity').slideDown();
        $('#examples_quiet').slideDown();
        $('#examples_precision').slideDown();
        $('#examples_solver').slideDown();
        
        $('h2.container_container_title').click(function(){     $(this).parent('div.body_content_container_container').find('div.content_content_box').slideToggle();
            //$(document.body).scrollTop($(this).offset().top);
        });

        $(document).find('div.body_content_container').find('div.content_box').slideUp();

        $('h1.container_title').click(function(){
            $(this).parent('div.body_content_container').find('div.content_box').slideToggle();
            //$(document.body).scrollTop($(this).offset().top);
        });
    }
    
    
});

function hide_button_press(){
    $("#hide_button_div").slideUp(0);
    hide_ig();
}
function apply_button_press(){
    $("#hide_button_div").slideDown(0);
    apply_options();
}

function hide_ig() {
    $("#ig_polyn_polynom_notcont_without").slideUp();
    $("#ig_polyn_polynom_notcont_with").slideUp();
    $("#ig_polyn_spline_notcont_without").slideUp();
    $("#ig_polyn_spline_cont_without").slideUp();
    $("#ig_polyn_spline_diff_without").slideUp();
    $("#ig_polyn_spline_cont_with").slideUp();
    $("#ig_polyn_spline_diff_with").slideUp();
    $("#ig_spline_polynom_notcont_without").slideUp();
    $("#ig_spline_polynom_notcont_with").slideUp();
    $("#ig_spline_spline_notcont_without").slideUp();
    $("#ig_spline_spline_notcont_with").slideUp();
    $("#ig_spline_spline_cont_without").slideUp();
    $("#ig_spline_spline_diff_without").slideUp();
    $("#ig_spline_spline_cont_with").slideUp();
    $("#ig_spline_spline_diff_with").slideUp();
    $("#ig_polyn_spline_notcont_with").slideUp();
}
function apply_options() {
    
    hide_ig();
    
    function_menu = $("#ig_function_menu");
    approximant_menu = $("#ig_approximant_menu");
    smoothness_menu = $("#ig_smoothness_menu");
    interpolation_menu = $("#ig_interpolation_menu");
    
    if (function_menu.val() == 1) {
        if (approximant_menu.val() == 1) {
            if (smoothness_menu.val() == 0) {
                if (interpolation_menu.val() == 0) {
                    //function(polynomial)
                    //approximant(polynomial)
                    //smoothness(continuous)
                    //interpolation(with)
                    $("#ig_polyn_polynom_notcont_with").slideDown();
                }
                else if (interpolation_menu.val() == 1){
                    //function(polynomial)
                    //approximant(polynomial)
                    //smoothness(continuous)
                    //interpolation(without)
                    $("#ig_polyn_polynom_notcont_without").slideDown();
                }
            }
            else if (smoothness_menu.val() == 1){
                if (interpolation_menu.val() == 0) {
                    //function(polynomial)
                    //approximant(polynomial)
                    //smoothness(differentiable)
                    //interpolation(with)
                    $("#ig_polyn_polynom_notcont_with").slideDown();
                }
                else {
                    //function(polynomial)
                    //approximant(polynomial)
                    //smoothness(differentiable)
                    //interpolation(without)
                    $("#ig_polyn_polynom_notcont_without").slideDown();
                }
            }
            else {
                if (interpolation_menu.val() == 0) {
                    //function(polynomial)
                    //approximant(polynomial)
                    //smoothness(non-continuous)
                    //interpolation(with)
                    $("#ig_polyn_polynom_notcont_with").slideDown();
                }
                else {
                    //function(polynomial)
                    //approximant(polynomial)
                    //smoothness(non-continuous)
                    //interpolation(without)
                    $("#ig_polyn_polynom_notcont_without").slideDown();
                }
            }
        }
        else if (approximant_menu.val() == 2){
            if (smoothness_menu.val() == 0) {
                if (interpolation_menu.val() == 0) {
                    //function(polynomial)
                    //approximant(spline)
                    //smoothness(continuous)
                    //interpolation(with)
                    $("#ig_polyn_spline_cont_with").slideDown();
                }
                else {
                    //function(polynomial)
                    //approximant(spline)
                    //smoothness(continuous)
                    //interpolation(without)
                    $("#ig_polyn_spline_cont_without").slideDown();
                }
            }
            else if (smoothness_menu.val() == 1){
                if (interpolation_menu.val() == 0) {
                    //function(polynomial)
                    //approximant(spline)
                    //smoothness(differentiable)
                    //interpolation(with)
                    $("#ig_polyn_spline_diff_with").slideDown();
                }
                else {
                    //function(polynomial)
                    //approximant(spline)
                    //smoothness(differentiable)
                    //interpolation(without)
                    $("#ig_polyn_spline_diff_without").slideDown();
                }
            }
            else {
                if (interpolation_menu.val() == 0) {
                    //function(polynomial)
                    //approximant(spline)
                    //smoothness(non-continuous)
                    //interpolation(with)
                    $("#ig_polyn_spline_notcont_with").slideDown();
                }
                else {
                    //function(polynomial)
                    //approximant(spline)
                    //smoothness(non-continuous)
                    //interpolation(without)
                    $("#ig_polyn_spline_notcont_without").slideDown();
                }
            }
        }
        else {
            //no approximant type selected
        }
    }
    else if (function_menu.val() == 2){
        if (approximant_menu.val() == 1) {
            if (smoothness_menu.val() == 0) {
                if (interpolation_menu.val() == 0) {
                    //function(spline)
                    //approximant(polynomial)
                    //smoothness(continuous)
                    //interpolation(with)
                    $("#ig_spline_polynom_notcont_with").slideDown();
                }
                else {
                    //function(spline)
                    //approximant(polynomial)
                    //smoothness(continuous)
                    //interpolation(without)
                    $("#ig_spline_polynom_notcont_without").slideDown();
                }
            }
            else if (smoothness_menu.val() == 1){
                if (interpolation_menu.val() == 0) {
                    //function(spline)
                    //approximant(polynomial)
                    //smoothness(differentiable)
                    //interpolation(with)
                    $("#ig_spline_polynom_notcont_with").slideDown();
                }
                else {
                    //function(spline)
                    //approximant(polynomial)
                    //smoothness(differentiable)
                    //interpolation(without)
                    $("#ig_spline_polynom_notcont_without").slideDown();
                }
            }
            else {
                if (interpolation_menu.val() == 0) {
                    //function(spline)
                    //approximant(polynomial)
                    //smoothness(non-continuous)
                    //interpolation(with)
                    $("#ig_spline_polynom_notcont_with").slideDown();
                }
                else {
                    //function(spline)
                    //approximant(polynomial)
                    //smoothness(non-continuous)
                    //interpolation(without)
                    $("#ig_spline_polynom_notcont_without").slideDown();
                }
            }
        }
        else if (approximant_menu.val() == 2){
            if (smoothness_menu.val() == 0) {
                if (interpolation_menu.val() == 0) {
                    //function(spline)
                    //approximant(spline)
                    //smoothness(continuous)
                    //interpolation(with)
                    $("#ig_spline_spline_cont_with").slideDown();
                }
                else {
                    //function(spline)
                    //approximant(spline)
                    //smoothness(continuous)
                    //interpolation(without)
                    $("#ig_spline_spline_cont_without").slideDown();
                }
            }
            else if (smoothness_menu.val() == 1){
                if (interpolation_menu.val() == 0) {
                    //function(spline)
                    //approximant(spline)
                    //smoothness(differentiable)
                    //interpolation(with)
                    $("#ig_spline_spline_diff_with").slideDown();
                }
                else {
                    //function(spline)
                    //approximant(spline)
                    //smoothness(differentiable)
                    //interpolation(without)
                    $("#ig_spline_spline_diff_without").slideDown();
                }
            }
            else {
                if (interpolation_menu.val() == 0) {
                    //function(spline)
                    //approximant(spline)
                    //smoothness(non-continuous)
                    //interpolation(with)
                    $("#ig_spline_spline_notcont_with").slideDown();
                }
                else {
                    //function(spline)
                    //approximant(spline)
                    //smoothness(non-continuous)
                    //interpolation(without)
                    $("#ig_spline_spline_notcont_without").slideDown();
                }
            }
        }
        else {
            //no approximant type selected
        }
    }
    else {
        //no function type selected
    }


};
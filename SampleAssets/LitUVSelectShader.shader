Shader "Custom/LitUVSelectShader"
{
    Properties
    {
        _MainTex("Texture", 2D) = "white" {}
        _UseUV1("Use UV1 (0 = UV0, 1 = UV1)", Float) = 0
    }
        SubShader
        {
            Tags { "RenderType" = "Opaque" }
            LOD 200

            CGPROGRAM
            #pragma surface surf Standard fullforwardshadows

            sampler2D _MainTex;
            float _UseUV1;

            struct Input
            {
                float2 uv_MainTex;        // Uses TEXCOORD0
                float2 uv2_MainTex;
            };

            void surf(Input IN, inout SurfaceOutputStandard o)
            {
                fixed4 c = tex2D(_MainTex, IN.uv_MainTex);
                fixed4 c1 = tex2D(_MainTex, IN.uv2_MainTex);
                fixed4 cf= lerp(c, c1, _UseUV1);
                o.Albedo = cf.rgb;
                o.Alpha = cf.a;
            }
            ENDCG
        }
            FallBack "Diffuse"
}

Shader "Custom/UVSelectShader"
{
    Properties
    {
        _MainTex("Texture", 2D) = "white" {}
        _UseUV1("Use UV1 (0 = UV0, 1 = UV1)", Float) = 0
    }
    SubShader
    {
        Tags { "RenderType" = "Opaque" }
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #include "UnityCG.cginc"

            sampler2D _MainTex;
            float _UseUV1;

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
                float2 uv1 : TEXCOORD1;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float2 uv1 : TEXCOORD1;
                float4 vertex : SV_POSITION;
            };

            v2f vert(appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;// lerp(v.uv, v.uv1, step(0.5, _UseUV1));
                o.uv1 = v.uv1;
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                float4 tex1Color= tex2D(_MainTex, i.uv);
                float4 tex2Color = tex2D(_MainTex, i.uv1);
                return lerp(tex1Color, tex2Color, _UseUV1);
            }
            ENDCG
        }
    }
}
